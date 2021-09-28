import xesmf as xe
import xarray as xr


def regrid_dataset(
    path_to_input_dataset, path_to_output_grid, in_config, out_config, path_to_output
):
    input_ds = xr.open_dataset(path_to_input_dataset).rename({
        in_config["lat-name"]: "lat", in_config["lon-name"]: "lon"
    })
    new_grid = xr.open_dataset(path_to_output_grid).rename({
        out_config["lat-name"]: "lat", out_config["lon-name"]: "lon"
    })

    # ASSUME: bilinear regridding (recommended by xesmf)
    regridder = xe.Regridder(input_ds, new_grid, "bilinear")

    regridded_data = regridder(input_ds)

    # Reattach non-lat/lon coordinates that have been lost
    for coord in [out_config["x-name"], out_config["y-name"]]:
        regridded_data.coords[coord] = new_grid.coords[coord]

    regridded_data.to_netcdf(path_to_output)


if __name__ == "__main__":
    regrid_dataset(
        path_to_input_dataset=snakemake.input.in_data,
        path_to_output_grid=snakemake.input.out_data,
        in_config=snakemake.params.in_config,
        out_config=snakemake.params.out_config,
        path_to_output=snakemake.output[0]
    )
