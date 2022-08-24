import os
import pandas as pd
import geopandas as gpd
import numpy as np

import sys
sys.path.append(f'{os.getcwd()}/src')
from utils import calculate_intersection_area, create_square_grids_from_centroids

print('Assigning municipalities to cells')
print('Note: the full script will take about 30-40 minutes to run')

processed_cells = gpd.read_file('/data/processed_data/cells/cells.geojson')
municipalities = gpd.read_file('/data/processed_data/municipalities/municipalities.geojson')

## Create grid polygons
# Create square grids from cell centroids
cells_grid = create_square_grids_from_centroids(processed_cells, cell_size=100, crs='epsg:25832')

## Spatial join
print('Start spatial join, which will take about 30 minutes')
# Spatial join grid cells with municipalties
cells_join_muni = gpd.sjoin(
    cells_grid, municipalities[['geometry','Kode']], how='left', predicate='intersects'
).drop(['index_right'], axis=1)
print('Finished spatial join')

## Handle duplicated cells
cells_join_muni_notnull = cells_join_muni[cells_join_muni.Kode.notnull()]
cells_join_muni_dup = cells_join_muni_notnull[cells_join_muni_notnull.ddkncelle100m.duplicated(keep=False)]
cells_join_muni_nondup = cells_join_muni_notnull[~cells_join_muni_notnull.ddkncelle100m.duplicated(keep=False)]
assert len(cells_join_muni_dup)+len(cells_join_muni_nondup) == len(cells_join_muni_notnull)

# For cells joined with more than one municipalities, assign it to the largest overlapped muni
# Get municipalities geometry
cells_join_muni_dup = cells_join_muni_dup.merge(
    municipalities[['Kode','geometry']].rename(columns={'geometry':'geometry_muni'}),
    how='left', on='Kode'
)

# Calculate intersected area between cell grids and muni polygons
cells_join_muni_dup['intersected_area'] = cells_join_muni_dup.apply(
    lambda row: calculate_intersection_area(row['geometry'], row['geometry_muni']), axis=1
)

# Sort by area so largest area is last, then drop duplicates and keep last/largest
cells_join_muni_dup_fixed = cells_join_muni_dup.sort_values(
    ['ddkncelle100m','intersected_area']
).drop_duplicates('ddkncelle100m', keep='last')
assert len(cells_join_muni_dup_fixed) + len(cells_join_muni_nondup) == cells_join_muni_notnull.ddkncelle100m.nunique()

cells_join_muni_final = pd.concat([
    cells_join_muni_dup_fixed[list(cells_join_muni_nondup.columns)],
    cells_join_muni_nondup
]).sort_values('ddkncelle100m').reset_index(drop=True)
assert len(cells_join_muni_final) == cells_join_muni_final.ddkncelle100m.nunique()
assert cells_join_muni_final.Kode.isnull().sum() == 0

cells_join_muni_final.Kode = cells_join_muni_final.Kode.astype(int)

cells_join_muni_final.drop(['geometry'], axis=1).to_csv(
    '/data/processed_data/cells/cells_with_municipalities.csv', index=False
)
