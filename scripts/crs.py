from pathlib import Path

import cartopy.crs as ccrs
import xarray as xr


def get_cosmo_rea_crs(pole_latitude, pole_longitude):
    return ccrs.RotatedPole(
        pole_latitude=pole_latitude, pole_longitude=pole_longitude
    ).proj4_init


def get_newa_crs(path_to_ds):
    """
    ASSUME: we don't have to use this function because we have crs attribute info in the downloaded dataset
    pyproj.Proj(
        proj='lcc', # projection type: Lambert Conformal Conic
        lat_1=ds.TRUELAT1, lat_2=ds.TRUELAT2, # Cone intersects with the sphere
        lat_0=ds.MOAD_CEN_LAT, lon_0=ds.STAND_LON, # Center point
        a=6370000, b=6370000
    ).crs
    """
    ds = xr.open_dataset(path_to_ds)
    return ds.crs.proj


if __name__ == "__main__":
    if snakemake.params.data_source.startswith("newa"):
        crs = get_newa_crs(snakemake.input.data)
    if snakemake.params.data_source.startswith("cosmo"):
        crs = get_cosmo_rea_crs(
            snakemake.params.pole_latitude, snakemake.params.pole_longitude
        )

    Path(snakemake.output[0]).write_text(crs)
