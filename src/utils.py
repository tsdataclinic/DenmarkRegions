import pandas as pd
import geopandas as gpd
import numpy as np
import libpysal
import spopt
from spopt.region import MaxPHeuristic as MaxP
from shapely.geometry import Polygon


def calculate_intersection_area(p1, p2):
    """
    Calculate the intersected areas between two polygons. The unit will be the same as what the
    input polygons are using.

    p1: shapely.geometry.Polygon
        polygon 1
    p2: shapely.geometry.Polygon
        polygon 2

    Returns
    -------
    GeoSeries.area
    """
    return p1.intersection(p2).area


def create_square_grids_from_centroids(
    gdf, cell_size=100, crs="EPSG:25832", crs_final="EPSG:4326"
):
    """
    Create square grid cells from centroid coordinates.

    gdf: gpd.GeoDataFrame
        geo dataframe of interest
    cell_size: int or float
        the length of one side of the cell, the unit should align with crs specified
    crs: str
        epsg string, cell_size should have the same unit as the crs specified
    crs_final: str
        epsg string, the dataframe returned will be projected to this crs

    Returns
    -------
    gpd.GeoDataFrame
    """
    gdf_meter = gdf.to_crs(crs)
    gdf_meter["centroid_x"] = gdf_meter.geometry.x
    gdf_meter["centroid_y"] = gdf_meter.geometry.y

    polygons = []
    f = int(cell_size / 2)
    for idx, row in gdf_meter.iterrows():
        x, y = row.centroid_x, row.centroid_y
        polygons.append(
            Polygon([(x - f, y + f), (x + f, y + f), (x + f, y - f), (x - f, y - f)])
        )

    # switch geometry to grid polygons
    gdf_meter_grid = gpd.GeoDataFrame(gdf_meter, geometry=polygons, crs=crs).drop(
        ["centroid_x", "centroid_y"], axis=1
    )

    # crs projection
    gdf_grid = (
        gdf_meter_grid.to_crs(crs_final) if crs != crs_final else gdf_meter_grid.copy()
    )

    return gdf_grid


def get_poly_with_pop(poly, cells):
    """
    Get population and household counts for each polygon by joinning cells data
    with polygons. When the centroid of the cell falls into the polygon, all cell
    values will be assign to that polygon.

    poly: gpd.GeoDataFrame
        geo dataframe of all the polygons
    cells: gpd.GeoDataFrame
        geo dataframe of all the cells with the population and household values

    Returns
    -------
    gpd.GeoDataFrame
    """
    poly = poly.to_crs("epsg:4326")
    poly["id"] = poly.index
    geoms = gpd.points_from_xy(cells.lng, cells.lat)
    geoCells = gpd.GeoDataFrame(cells, geometry=geoms, crs="EPSG:4326")
    cells_w_poly = gpd.sjoin(geoCells, poly, how="left")
    poly_stats_all = cells_w_poly.groupby("id").sum().reset_index()

    hust_pers_cols = poly_stats_all.filter(regex=r"(hust|pers)").columns.tolist()
    polygon_w_pop = poly.merge(
        poly_stats_all[["id"] + hust_pers_cols], how="left", on="id"
    )

    hust_cols = [x for x in polygon_w_pop.columns if "hust" in x]
    pers_cols = [x for x in polygon_w_pop.columns if "pers" in x]
    polygon_w_pop["min_hust"] = polygon_w_pop[hust_cols].apply(min, axis=1)
    polygon_w_pop["min_pers"] = polygon_w_pop[pers_cols].apply(min, axis=1)

    polygon_w_pop["min_val"] = [
        min(x, y / 2) for x, y in zip(polygon_w_pop.min_hust, polygon_w_pop.min_pers)
    ]
    polygon_w_pop["min_2020"] = [
        min(x, y / 2) for x, y in zip(polygon_w_pop.hust2020, polygon_w_pop.pers2020)
    ]
    polygon_w_pop[["min_val", "min_2020"]] = polygon_w_pop[
        ["min_val", "min_2020"]
    ].fillna(0)

    return polygon_w_pop


def get_poly_with_interpolated_pop(
    poly, cells, cell_size=100, original_crs="EPSG:25832", new_crs="EPSG:4326"
):
    """
    Get population and household counts for each polygon by joinning cells data
    with polygons. When a cell overlaps with more than one polygon, the intersected
    areas and an overlap ratio will be calculated, and the cell values will be first
    multiplied by the ratio and then added to the polygon (areal interpolation).

    poly: gpd.GeoDataFrame
        geo dataframe of all the polygons
    cells: gpd.GeoDataFrame
        geo dataframe of all the cells with the population and household values
    cell_size: int or float
        the length of one side of the cell, the unit should align with crs specified
    original_crs: str
        epsg string, cell_size should have the same unit as the crs specified
    new_crs: str
        epsg string, the dataframe returned will be projected to this crs

    Returns
    -------
    gpd.GeoDataFrame
    """
    poly = poly.to_crs(new_crs)
    poly["id"] = poly.index

    cell_grid = create_square_grids_from_centroids(
        cells, cell_size=cell_size, crs=original_crs, crs_final=new_crs
    )
    cells_join_poly = gpd.sjoin(poly, cell_grid, how="left").drop(
        ["index_right"], axis=1
    )

    # find intersected areas of cells and polygons
    cell_grid["geometry_cell"] = cell_grid["geometry"]
    cells_join_poly_interpolate = cells_join_poly[
        cells_join_poly.ddkncelle100m.notnull()
    ].merge(
        cell_grid[["ddkncelle100m", "geometry_cell"]], how="left", on="ddkncelle100m"
    )

    cells_join_poly_interpolate["intersected_area"] = cells_join_poly_interpolate.apply(
        lambda row: calculate_intersection_area(row["geometry"], row["geometry_cell"]),
        axis=1,
    )

    # grid ratio: intersected_area/grid area
    cells_join_poly_interpolate["grid_ratio"] = (
        cells_join_poly_interpolate.intersected_area
        / cells_join_poly_interpolate.geometry_cell.area
    )

    # multiply household and pop counts with grid ratio
    hust_pers_cols = cells_join_poly_interpolate.filter(
        regex=r"(hust|pers)"
    ).columns.tolist()
    for c in hust_pers_cols:
        cells_join_poly_interpolate[f"{c}_interpolated"] = (
            cells_join_poly_interpolate[c] * cells_join_poly_interpolate["grid_ratio"]
        )

    # add interpolated poly stats
    interpolated_cols = [c for c in cells_join_poly_interpolate if "interpolated" in c]
    poly_stats = (
        cells_join_poly_interpolate.groupby("id")[interpolated_cols]
        .sum(min_count=1)
        .reset_index()
    )
    poly_w_pop = poly.merge(poly_stats, how="left", on="id").fillna(0)

    # get min values
    hust_cols = [x for x in poly_w_pop.columns if "hust" in x]
    pers_cols = [x for x in poly_w_pop.columns if "pers" in x]
    poly_w_pop["min_hust"] = poly_w_pop[hust_cols].apply(min, axis=1)
    poly_w_pop["min_pers"] = poly_w_pop[pers_cols].apply(min, axis=1)

    poly_w_pop["min_val"] = [
        min(x, y / 2) for x, y in zip(poly_w_pop.min_hust, poly_w_pop.min_pers)
    ]
    poly_w_pop["min_2020"] = [
        min(x, y / 2)
        for x, y in zip(
            poly_w_pop.hust2020_interpolated, poly_w_pop.pers2020_interpolated
        )
    ]
    poly_w_pop[["min_val", "min_2020"]] = poly_w_pop[["min_val", "min_2020"]].fillna(0)

    return poly_w_pop


def get_buffered_weights(poly_w_pop, weight="Rook"):
    """
    Helper function to get buffered weights for each component.

    poly_w_pop: gpd.GeoDataFrame
        geo dataframe of the polygons with household and population counts
    weight: str
        name of the contiguity-based spatial weight, such as Queen or Rook

    Returns
    -------
    poly_w_pop: gpd.GeoDataFrame
    w_new: pysal spatial weights class
    """

    if weight == "Rook":
        w = libpysal.weights.Rook
    elif weight == "Queen":
        w = libpysal.weights.Queen
    else:
        exc = "Only support Rook or Queen weight now"
        raise Exception(exc)

    w_a = w.from_dataframe(poly_w_pop)
    poly_w_pop["component"] = w_a.component_labels
    component_sums = poly_w_pop.groupby("component").sum().reset_index()
    component_counts = (
        poly_w_pop.groupby("component").count()[["id"]] * 100 / poly_w_pop.shape[0]
    )
    component_counts = component_counts.reset_index()
    empty_comps = component_sums[component_sums.min_pers == 0].component.values

    # remove empty comps
    poly_w_pop = poly_w_pop[~poly_w_pop.component.isin(empty_comps)]
    poly_w_pop = poly_w_pop.reset_index()
    poly_w_pop["id"] = poly_w_pop.index
    w_a = w.from_dataframe(poly_w_pop)

    # buffer small components
    small_components = component_counts[component_counts.id < 25].component.values
    poly_to_buffer = poly_w_pop[poly_w_pop.component.isin(small_components)]
    poly_old = poly_w_pop[~poly_w_pop.component.isin(small_components)]
    w_b = libpysal.weights.fuzzy_contiguity(poly_w_pop, buffering=True, buffer=0.00015)

    unbuffered_neighbors = {k: w_a.neighbors[k] for k in poly_old.id.values}
    buffered_neighbors = {k: w_b.neighbors[k] for k in poly_to_buffer.id.values}
    new_neighbors = unbuffered_neighbors
    new_neighbors.update(buffered_neighbors)

    unbuffered_weights = {k: w_a.weights[k] for k in poly_old.id.values}
    buffered_weights = {k: w_b.weights[k] for k in poly_to_buffer.id.values}
    new_weights = unbuffered_weights
    new_weights.update(buffered_weights)

    w_new = libpysal.weights.W(new_neighbors, new_weights)
    poly_w_pop["component"] = w_new.component_labels
    new_poly_stats = poly_w_pop.groupby("component").sum().reset_index()
    comps_to_remove = new_poly_stats[new_poly_stats.min_val < 50].component.values
    if len(comps_to_remove) > 0:
        poly_w_pop = poly_w_pop[~poly_w_pop.component.isin(comps_to_remove)]
        poly_w_pop, w_new = get_buffered_weights(poly_w_pop)

    return poly_w_pop, w_new
