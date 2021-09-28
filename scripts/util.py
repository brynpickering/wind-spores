import numpy as np


def coeffs_to_wind_speed(coeff_ds, height):
    return coeff_ds.A * (np.log(height) - coeff_ds.log_z)


def reassign_lat_lon_names(input_ds, output_ds, lat_name=None, lon_name=None):
    for var in [lon_name, lat_name]:
        if var is not None:
            output_ds.coords[var] = input_ds.coords[var]
    return output_ds
