import xarray as xr

import util


def coeff_to_ws(path_to_coeffs, height, lon_name, lat_name, path_to_output):
    coeffs = xr.open_dataset(path_to_coeffs)

    wind_speed = util.coeffs_to_wind_speed(coeffs, float(height))

    util.reassign_lat_lon_names(coeffs, wind_speed, lat_name, lon_name)

    wind_speed.to_netcdf(path_to_output)


if __name__ == "__main__":
    coeff_to_ws(
        path_to_coeffs=snakemake.input.coeffs,
        height=snakemake.wildcards.height,
        lon_name=snakemake.params.lon_name,
        lat_name=snakemake.params.lat_name,
        path_to_output=snakemake.output[0]
    )
