import xarray as xr
import pandas as pd


def slice_cosmo_dataset(path_to_input, path_to_cosmo_coords, bounding_box, cosmo_config, year, path_to_output):
    input_ds = xr.open_dataset(path_to_input)
    coords = xr.open_dataset(path_to_cosmo_coords)
    input_ds.coords[cosmo_config["lat-name"]] = coords.coords[cosmo_config["lat-name"]]
    input_ds.coords[cosmo_config["lon-name"]] = coords.coords[cosmo_config["lon-name"]]

    bounded_ds = input_ds.where(
        (input_ds[cosmo_config["lat-name"]] >= bounding_box["lat-range"][0]) &
        (input_ds[cosmo_config["lat-name"]] <= bounding_box["lat-range"][1]) &
        (input_ds[cosmo_config["lon-name"]] >= bounding_box["lon-range"][0]) &
        (input_ds[cosmo_config["lon-name"]] <= bounding_box["lon-range"][1]),
        drop=True
    )

    # Align the timestamp of a particular observation with our preference: start of hour
    # (input is based on end of observation)
    bounded_ds.coords["time"] = bounded_ds.time.to_index() + pd.Timedelta(-1, unit="h")
    bounded_ds.loc[{"time": str(year)}].to_netcdf(path_to_output)


if __name__ == "__main__":
    slice_cosmo_dataset(
        path_to_input=snakemake.input.wind_speed,
        path_to_cosmo_coords=snakemake.input.cosmo_coords,
        bounding_box=snakemake.params.bounding_box,
        cosmo_config=snakemake.params.cosmo_config,
        year=snakemake.wildcards.year,
        path_to_output=snakemake.output[0]
    )
