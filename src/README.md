# Code Structure

This folder contains the main ETL scripts and notebooks that process the data and run the clustering analysis. Below are descriptions of the key code files.

## ETL scripts

- `cells`
    - `process_cells_data.py`: Process raw cells csv file, including imputing missing population and household values, and shifting cells coordinates to represent the centroids
    - `assign_municipalities_to_cells.py`: Spatial join processed cells with municipalities. For cells intersected with more than one municipality, assign it to the largest overlapped municipality
- `municipalities`
    - `process_municipalities.py`: Process raw municipalities file
- `cluster`
    - `cluster_polygons.py`: Run Max-P clustering on polygons. Population and household values are not interpolated
    - `cluster_polygons_w_areal_interpolation.py`: Run Max-P clustering on polygons. Population and household values will be interpolated
    - `cluster_gridded_polygons.py`: Run Max-P clustering on gridded polygons. 


## Notebooks
- `notebooks` 
    - `calculate_metrics.ipynb`: Contains tables and graphs comparing key metrics of different clustering results
    - `Example.ipynb`: Contains an example of the clustering process using sample data available in `data/`