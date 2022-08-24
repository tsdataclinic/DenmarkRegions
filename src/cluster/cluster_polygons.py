import os
import pandas as pd
import numpy as np
import geopandas as gpd
import spopt
from spopt.region import MaxPHeuristic as MaxP
import matplotlib.pyplot as plt

import sys
sys.path.append(f'{os.getcwd()}/src')
from utils import get_poly_with_pop, get_buffered_weights

results_dir = '/data/results/polygons_compact_rook/'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
    os.makedirs(results_dir+'plots/')
    os.makedirs(results_dir+'geometries/')

cells = pd.read_csv("/data/processed_data/cells/cells_with_municipalities.csv")
geoCells = gpd.GeoDataFrame(cells, geometry= gpd.points_from_xy(cells.x_bin,cells.y_bin), crs='EPSG:25832')
geoCells = geoCells.to_crs('epsg:4326')
cells_to_cluster = geoCells.assign(
    lat = geoCells.geometry.y, 
    lng = geoCells.geometry.x,
)
municipalities = cells[~cells.Kode.isna()].Kode.unique()
municipalities = [int(m) for m in municipalities]

clus_res_dfs = []
municipalities_w_issues = []

for m in municipalities:
    print(m)
    try:
        previous_dir = '/data/Polygons/polygons_cleaned'
        current_dir = '/data/data_org_polygons/data/processed/polygons/polygons_compact/polygon_compact'
        poly = gpd.read_file(f'{current_dir}_{m}.geojson')
        if 'munic_name' in poly.columns:
            poly = poly.drop(columns=['munic_name'])

        poly = poly.to_crs('epsg:4326')

        to_cluster = cells_to_cluster[cells_to_cluster.Kode==m]

        poly_w_pop = get_poly_with_pop(poly, to_cluster)

        poly_to_cluster, w = get_buffered_weights(poly_w_pop)

        threshold = 50
        top_n = 2
        print('constructing model')
        model = MaxP(poly_to_cluster, w, [], 'min_val', threshold, top_n)
        model.solve()

        poly_to_cluster['max_p'] = model.labels_
        clus_res_dfs.append(poly_to_cluster[['id','min_hust','min_pers','min_val','max_p','geometry']].assign(Kode=m))

        regions = poly_to_cluster.drop(['index','id','munic_code','component'],axis=1).dissolve(by='max_p', aggfunc='sum')
        regions.to_file(results_dir+'geometries/Kode_{}.geojson'.format(m),driver='GeoJSON')

        fig,ax = plt.subplots(figsize=(20,15))
        p = poly_to_cluster.plot(column='max_p',ax=ax, cmap="Spectral", edgecolor='white')
        ax.title.set_text('Kode_{}'.format(m))
        plt.savefig(results_dir+'plots/Kode_{}.png'.format(m))
    except Exception as e:
        print(e)
        municipalities_w_issues.append(m)

cluster_results = pd.concat(clus_res_dfs)

cluster_results.drop(columns='geometry').to_csv(results_dir+'polygon_cluster_mapping.csv',index=False)
cluster_results.to_file(results_dir+'polygon_cluster_mapping.geojson',driver='GeoJSON')

print("Municipalities with issues:", municipalities_w_issues)