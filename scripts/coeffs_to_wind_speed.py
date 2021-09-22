import xarray as xr

import util


def coeff_to_ws(path_to_coeffs, height, path_to_output):
    wind_speed = util.coeffs_to_wind_speed(xr.open_dataset(path_to_coeffs), float(height))
    wind_speed.to_netcdf(path_to_output)


if __name__ == "__main__":
    coeff_to_ws(
        path_to_coeffs=snakemake.input.coeffs,
        height=snakemake.wildcards.height,
        path_to_output=snakemake.output[0]
    )
