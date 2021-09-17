import xarray as xr
import numpy as np
import pandas as pd

import util


def speed_to_capacityfactor(
    path_to_wind_speed_coeffs, path_to_power_curves, turbine_name, turbine_config, path_to_output
):
    wind_speed_coeffs = xr.open_dataset(path_to_wind_speed_coeffs)
    wind_speed = util.coeffs_to_wind_speed(wind_speed_coeffs, turbine_config["height"])

    power_curve = pd.read_csv(path_to_power_curves, index_col=0).loc[:, turbine_config["long-name"]]

    power = xr.apply_ufunc(
        interpolate_along_power_curve, wind_speed, kwargs={"power_curve": power_curve}
    )

    cf_ds = xr.Dataset({
        "capacityfactor": power.expand_dims({"turbine": pd.Index([turbine_name])}).astype(np.float32)
    })

    cf_ds.to_netcdf(path_to_output)


def interpolate_along_power_curve(array, power_curve):
    return np.interp(array, power_curve.index, power_curve.values)


if __name__ == "__main__":
    speed_to_capacityfactor(
        path_to_wind_speed_coeffs=snakemake.input.wind_speed_coeffs,
        path_to_power_curves=snakemake.input.power_curves,
        turbine_name=snakemake.wildcards.turbine_name,
        turbine_config=snakemake.params.turbine_config,
        path_to_output=snakemake.output[0]
    )
