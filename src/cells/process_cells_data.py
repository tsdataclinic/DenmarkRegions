import pandas as pd 
import geopandas as gpd
import numpy as np

raw_data_dir = '/data'
processed_data_dir = '/data/processed_data'

## Load data
cells = pd.read_csv(f'{raw_data_dir}/Hectare_cells/grid_distribution_v3.txt')

# 99999 means the value could be 1 or 2, to be imputed later
cells = cells.replace(99999, np.nan)

## Imputation
# transform the cells df to two long format tables and then combine them
pers = pd.melt(cells, id_vars=['ddkncelle100m'], value_vars=[c for c in cells if 'pers' in c])
pers['year'] = pers.variable.str.extract('(\d+)').astype(int)

hust = pd.melt(cells, id_vars=['ddkncelle100m'], value_vars=[c for c in cells if 'hust' in c])
hust['year'] = hust.variable.str.extract('(\d+)').astype(int)

cells_long = pd.merge(
    pers.rename(columns={'variable': 'pers_year', 'value': 'pers_original'}),
    hust.rename(columns={'variable': 'hust_year', 'value': 'hust_original'}),
    how='left', on=['ddkncelle100m', 'year']
)

# to be conservative: fill all missing person counts with 1
cells_long['pers'] = cells_long['pers_original'].fillna(1)

# # optional way - if number of people is above a threshold, set household counts to be 2, set the rest to be 1
# threshold = 10
# cells_long['hust'] = np.where(
#     (cells_long.hust_original.isnull())&(cells_long.pers_original>threshold),
#     2,
#     cells_long.hust_original
# )
# cells_long['hust'] = cells_long['hust'].fillna(1)

# to be conservative: fill all missing household counts with 1
cells_long['hust'] = cells_long['hust_original'].fillna(1)

# verify
assert cells_long.hust.isnull().sum() == 0
assert cells_long.pers.isnull().sum() == 0

# transform back to wide format
pers_wide = cells_long.pivot(index='ddkncelle100m', columns='pers_year', values='pers').reset_index().rename_axis(None, axis=1)
hust_wide = cells_long.pivot(index='ddkncelle100m', columns='hust_year', values='hust').reset_index().rename_axis(None, axis=1)
cells_interpolated = pd.merge(pers_wide, hust_wide, how='left', on='ddkncelle100m')[list(cells.columns)]
assert len(cells_interpolated) == len(cells)

## Get geometry
# add bins
x_coords  = cells_interpolated.ddkncelle100m.str.split('_').str[2].astype(float)
y_coords  = cells_interpolated.ddkncelle100m.str.split('_').str[1].astype(float)

# the coordinates represent lower left corners of each cell, so shift it to the centroid
cell_size = 100
x_coords, y_coords = x_coords*100+cell_size/2, y_coords*100+cell_size/2
cells_interpolated = cells_interpolated.assign(x_bin = x_coords, y_bin=y_coords)

## cache results
geoCells = gpd.GeoDataFrame(cells_interpolated, geometry= gpd.points_from_xy(x_coords,y_coords), crs='EPSG:25832')
geoCells = geoCells.to_crs('epsg:4326')
geoCells.to_file(f'{processed_data_dir}/cells/cells.geojson', driver='GeoJSON')

cells_width_edges = geoCells.assign(
    lat = geoCells.geometry.y, 
    lng = geoCells.geometry.x,
    cell_id = geoCells.index
).drop(['geometry'],axis=1)
cells_width_edges.to_csv(f'{processed_data_dir}/cells/cells.csv', index=False)
