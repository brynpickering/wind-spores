import numpy as np


def coeffs_to_wind_speed(coeff_ds, height):
    return coeff_ds.A * (np.log(height) - coeff_ds.log_z)