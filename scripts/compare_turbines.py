import xarray as xr


def compare_turbines(
    paths_to_turbine_datasets, path_to_output
):

    all_turbine_da = xr.concat(
        [xr.open_dataarray(file) for file in paths_to_turbine_datasets], dim="turbine"
    )

    best_turbines = all_turbine_da.idxmax("turbine")

    best_turbines.rename("best_turbines").to_netcdf(path_to_output)


if __name__ == "__main__":
    compare_turbines(
        paths_to_turbine_datasets=snakemake.input.turbine_metrics,
        path_to_output=snakemake.output[0]
    )
