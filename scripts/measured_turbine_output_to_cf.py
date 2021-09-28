import xarray as xr
import pandas as pd
import geopandas as gpd


def turbine_output_to_cf(
    path_to_turbine_output, path_to_turbines_per_site, path_to_current_turbine_sites,
    turbine_config, site_id, turbine_id, path_to_output
):
    turbine_name = turbine_config["long-name"]
    turbine_site_data = (
        gpd.read_file(path_to_current_turbine_sites)
        .set_index(["site_id", "turbine"])
        .loc[(site_id, turbine_name)]
    )
    n_turbines = (
        pd.read_csv(path_to_turbines_per_site, index_col=[0, 1])
        .loc[(site_id, turbine_name)]
        .astype(float)
    )
    n_turbines.index = n_turbines.index.astype(int)

    site_capacity = n_turbines.mul(turbine_site_data.capacity)

    turbine_output = (
        xr.open_dataset(path_to_turbine_output)
        [site_id]
        .sel(turbine=turbine_id)
        .dropna("time")
    )

    turbine_cf = (
        turbine_output
        .groupby("time.year").apply(
            lambda x: x / site_capacity.reindex(x.time.dt.year)
        )
        .reset_coords(drop=True)
        .expand_dims({"dataset": ["Measured"], "site": [site_id], "turbine": [turbine_id]})
        .to_dataset(name="capacityfactor")
    )
    turbine_cf.to_netcdf(path_to_output)


if __name__ == "__main__":
    turbine_output_to_cf(
        path_to_turbine_output=snakemake.input.turbine_output,
        path_to_turbines_per_site=snakemake.input.turbines_per_site,
        path_to_current_turbine_sites=snakemake.input.current_turbine_sites,
        turbine_config=snakemake.params.turbine_config,
        site_id=snakemake.wildcards.turbine_site,
        turbine_id=snakemake.wildcards.turbine_id,
        path_to_output=snakemake.output[0]
    )
