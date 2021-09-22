import cartopy.crs as ccrs
import xarray as xr
import pandas as pd
import geopandas as gpd
import shapely.geometry


def points_to_polys(path_to_input, data_source, config, output_driver, path_to_output):
    ds = xr.open_dataset(path_to_input)
    points = pd.MultiIndex.from_product(
        [ds.coords[config["x-name"]].to_index(), ds.coords[config["y-name"]].to_index()],
    )
    geoms = [
        shapely.geometry.Point(point) for point in points
    ]
    if data_source.startswith("cosmo"):
        crs = get_cosmo_rea_crs(
            config["pole-latitude"], config["pole-longitude"]
        )
    elif data_source == "newa":
        crs = get_wrf_crs(ds)

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


def get_cosmo_rea_crs(pole_latitude, pole_longitude):
    return ccrs.RotatedPole(
        pole_latitude=pole_latitude, pole_longitude=pole_longitude
    ).proj4_init


def get_wrf_crs(ds):
    """
    ASSUME: we don't have to use this function because we have crs attribute info in the downloaded dataset
    pyproj.Proj(
        proj='lcc', # projection type: Lambert Conformal Conic
        lat_1=ds.TRUELAT1, lat_2=ds.TRUELAT2, # Cone intersects with the sphere
        lat_0=ds.MOAD_CEN_LAT, lon_0=ds.STAND_LON, # Center point
        a=6370000, b=6370000
    ).crs
    """
    return ds.crs.proj


if __name__ == "__main__":
    points_to_polys(
        path_to_input=snakemake.input.data,
        data_source=snakemake.params.data_source,
        config=snakemake.params.config,
        output_driver=snakemake.params.output_driver,
        path_to_output=snakemake.output[0]
    )
