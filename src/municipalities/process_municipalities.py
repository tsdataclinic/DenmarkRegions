import geopandas as gpd

municips = gpd.read_file('/data/Municipality_shape/current/munic_final.geojson')
municips = municips.to_crs('EPSG:4326')
municips['Kode'] = municips['munic_code']

municips.to_file('/data/processed_data/municipalities/municipalities.geojson', driver='GeoJSON')
