from pathlib import Path

import xarray as xr
import pandas as pd
import geopandas as gpd
import shapely.geometry


def points_to_polys(path_to_input, path_to_crs, config, output_driver, path_to_output):
    ds = xr.open_dataset(path_to_input)
    points = pd.MultiIndex.from_product(
        [ds.coords[config["x-name"]].to_index(), ds.coords[config["y-name"]].to_index()],
    )
    geoms = [
        shapely.geometry.Point(point) for point in points
    ]

    crs = Path(path_to_crs).read_text()
    geopoints = gpd.GeoSeries(geoms, index=points, crs=crs)

    geopolys = geopoints.buffer(get_gridsize(points) / 2).envelope

    # The safest is to store geoinformation as WGS84/EPSG4326.
    # e.g. GeoJSON won't store information about a non-standard CRS
    geopolys_epsg4326 = geopolys.to_crs("EPSG:4326")

    geopolys_epsg4326.to_file(path_to_output, driver=output_driver)


def get_gridsize(points):
    diff = points.to_frame().diff()
    most_likely_diff = diff.where(diff != 0).median()
    assert round(most_likely_diff.std(), 1) == 0.0

    return most_likely_diff.mean()


if __name__ == "__main__":
    points_to_polys(
        path_to_input=snakemake.input.data,
        path_to_crs=snakemake.input.crs,
        config=snakemake.params.config,
        output_driver=snakemake.params.output_driver,
        path_to_output=snakemake.output[0]
    )
