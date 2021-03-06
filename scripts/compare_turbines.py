import xarray as xr


def compare_turbines(
    paths_to_turbine_datasets, path_to_output
):

    all_turbine_da = xr.concat(
        [xr.open_dataarray(file) for file in paths_to_turbine_datasets], dim="turbine"
    )

    best_turbines = all_turbine_da.idxmax("turbine").to_dataset(name="turbine_max")
    best_turbines["turbine_min"] = all_turbine_da.idxmin("turbine")

    interseasonal_da = (
        abs(all_turbine_da.sel(period="winter") - all_turbine_da.sel(period="summer")) /
        all_turbine_da.sel(period="all")
    )

    best_turbines["interseasonal_turbine_max"] = interseasonal_da.idxmax("turbine")
    best_turbines["interseasonal_turbine_min"] = interseasonal_da.idxmin("turbine")

    best_turbines.to_netcdf(path_to_output)


if __name__ == "__main__":
    compare_turbines(
        paths_to_turbine_datasets=snakemake.input.turbine_metrics,
        path_to_output=snakemake.output[0]
    )
