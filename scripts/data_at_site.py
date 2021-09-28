from pathlib import Path

import xarray as xr
import geopandas as gpd


def data_at_site(
    path_to_crs, path_to_modelled_data, path_to_current_turbine_sites,
    x_name, y_name, site_id, dataset_name, path_to_output
):
    modelled_ds = xr.open_dataarray(path_to_modelled_data)

    crs = Path(path_to_crs).read_text()
    turbine_sites = (
        gpd.read_file(path_to_current_turbine_sites).dissolve("site_id")
    )
    site_reprojected = turbine_sites.to_crs(crs).loc[site_id]
    x_coord = site_reprojected.geometry.centroid.x
    y_coord = site_reprojected.geometry.centroid.y

    modelled_site_data = (
        interp_data_to_site(x_coord, y_coord, x_name, y_name, modelled_ds)
        .expand_dims({"dataset": [dataset_name], "site": [site_id]})
    )
    modelled_site_data.to_netcdf(path_to_output)


def interp_data_to_site(x_coord, y_coord, x_name, y_name, modelled_data):
    return (
        modelled_data
        .interp(**{x_name: x_coord, y_name: y_coord})
        .reset_coords(drop=True)
    )


if __name__ == "__main__":
    data_at_site(
        path_to_crs=snakemake.input.crs,
        path_to_modelled_data=snakemake.input.modelled_data,
        path_to_current_turbine_sites=snakemake.input.current_turbine_sites,
        x_name=snakemake.params.x_name,
        y_name=snakemake.params.y_name,
        site_id=snakemake.wildcards.turbine_site,
        dataset_name=snakemake.params.dataset_name,
        path_to_output=snakemake.output[0]
    )
