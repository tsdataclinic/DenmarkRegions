import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt 
import libpysal
import spopt
import os
from shapely.ops import cascaded_union
from spopt.region import MaxPHeuristic as MaxP
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

## Loading raw data
cells = pd.read_csv("/data/processed_data/cells/cells_with_municipalities.csv")
# polygons = gpd.read_file('/data/Polygons/polygons_v2.geojson')
municipalities_shapes = gpd.read_file('/data/processed_data/municipalities/municipalities.geojson')
grid_shapes = gpd.read_file('../notebooks/grid_final.geojson')

## Setting path for results
results_dir = '/data/results/combined_grid/'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
    os.makedirs(results_dir+'plots/')
    os.makedirs(results_dir+'geometries/')

    
municipalities = municipalities_shapes[~municipalities_shapes.Kode.isna()].Kode.unique()
municipalities = [int(m) for m in municipalities]

clus_res_dfs = []
compactness_res_dfs = []
municipalities_w_issues = []

for m in municipalities:
    print(m)
    try:
        ## Combining Grid Cells with polygons
        print('combining grid cells')
        grid_mapping = pd.read_csv('/data/processed_data/grid_mapping/grid_clusters_{}.csv'.format(m))
        grid_in_poly = grid_mapping.merge(grid_shapes[['cell','geometry']], how='left', on='cell')
        grid_in_poly = grid_in_poly.merge(cells,how='left',left_on='cell',right_on='ddkncelle100m')

        grid_in_poly['component'] = 0
        ids = list(grid_in_poly.cluster.unique())

        for i in ids:
            tmp = grid_in_poly[grid_in_poly.cluster==i]
            w = libpysal.weights.Rook.from_dataframe(tmp)
            grid_in_poly.loc[grid_in_poly.cluster==i,'component'] = w.component_labels

        grid_in_poly['id'] = grid_in_poly.cluster + grid_in_poly.component/100
        grid_in_poly[['cell','Kode','id']].to_csv(
                '/data/processed_data/polygons/combined_grid_cells/grid_id_mapping_{}.csv'.format(m),index=False)
        grid_in_poly = gpd.GeoDataFrame(grid_in_poly,geometry=grid_in_poly.geometry)
        grid_groups = grid_in_poly.dissolve(by='id',aggfunc='sum').reset_index()

        ## Calculating parameter for Max-P threshold
        hust_cols = [x for x in grid_groups.columns if 'hust' in x]
        pers_cols = [x for x in grid_groups.columns if 'pers' in x]

        grid_groups['min_hust'] = grid_groups[hust_cols].apply(min,axis=1)
        grid_groups['min_pers'] = grid_groups[pers_cols].apply(min,axis=1)

        grid_groups['min_val'] = [min(x,y/2) for x,y in zip(grid_groups.min_hust, grid_groups.min_pers)]
        grid_groups['min_2020'] = [min(x,y/2) for x,y in zip(grid_groups.hust2020, grid_groups.pers2020)]
        grid_groups[['min_val','min_2020']] = grid_groups[['min_val','min_2020']].fillna(0)

        grid_groups.to_file('/data/processed_data/polygons/combined_grid_cells/combined_grid_{}.geojson'.format(m),driver='GeoJSON')

        ## Initializing the adjecancy for Max-P
        w = libpysal.weights.Rook.from_dataframe(grid_groups)
        grid_groups['poly_comp'] = w.component_labels
        comp_vals = grid_groups.groupby('poly_comp').sum()[['min_val']].reset_index()
        components_to_cluster = list(comp_vals[comp_vals.min_val>=50].poly_comp)
        grid_groups_minus_islands = grid_groups[grid_groups.poly_comp.isin(components_to_cluster)]
        w_rook = libpysal.weights.Rook.from_dataframe(grid_groups_minus_islands)

        ## Solving Max-P regionalization
        threshold = 50
        top_n = 2
        print('constructing model')
        model = MaxP(grid_groups_minus_islands, w_rook, [], 'min_val', threshold, top_n)
        model.solve()

        ## Compiling results
        grid_groups_minus_islands['max_p'] = model.labels_

        clus_res_dfs.append(grid_groups_minus_islands[['id','min_hust','min_pers','min_val','max_p','geometry']])

        regions = grid_groups_minus_islands.dissolve(by='max_p', aggfunc='sum').reset_index()
        regions.to_file(results_dir+'geometries/Kode_{}.geojson'.format(m),driver='GeoJSON')

        ipqs = regions.area * 4 * np.pi / (regions.boundary.length**2)
        ipqs = ipqs.to_frame().reset_index().rename(columns={0:'compactness','max_p':'clus'})
        ipqs['method'] = 'combined_grid'
        ipqs['Kode'] = m
        compactness_res_dfs.append(ipqs)

        fig,ax = plt.subplots(figsize=(20,15))
        p = regions.plot(column='max_p',ax=ax, cmap="Spectral", edgecolor='white')
        ax.title.set_text('Kode_{}'.format(m))
        plt.savefig(results_dir+'plots/Kode_{}.png'.format(m))
    except Exception as e:
        print(e)
        municipalities_w_issues.append(m)

## Writing out results
cluster_results = pd.concat(clus_res_dfs)
compactness_results = pd.concat(compactness_res_dfs)

cluster_results.drop(columns='geometry').to_csv(results_dir+'combined_gridcell_cluster_mapping.csv',index=False)
compactness_results.to_csv(results_dir+'cluster_compactness.csv',index=False)
cluster_results.to_file(results_dir+'combined_gridcell_cluster_mapping.geojson',driver='GeoJSON')