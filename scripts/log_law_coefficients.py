from sklearn.linear_model import LinearRegression
import numpy as np
import pandas as pd
import xarray as xr


def wind_speed_to_coeffs(path_to_input, height_weights, x_name, y_name, path_to_output):
    height_weights_df = pd.Series(height_weights)
    A_z_ds = get_A_z(path_to_input, x_name, y_name, height_weights_df)

    A_z_ds.to_netcdf(path_to_output)


def get_A_z(path_to_input, x_name, y_name, height_weights):
    """
    Apply log law regression to data for a specific year.

    Parameters
    ----------
    path_to_input : str
        points to a single NEWA wind file.
    x_name, y_name : str
        names for the geographic coordinates used
        (may or may not be reprojected from original data)
    height_weights : pandas.DataFrame
        index = height, columns = weight corresponding to each height.

    Returns
    -------
    A and z log-law coefficients.

    """
    # Get list of subdirectories

    wind_speeds = xr.open_dataset(path_to_input).WS

    wind_speed_stacked = wind_speeds.stack(x_y_time=[x_name, y_name, "time"])

    A, log_z = extrapolate(wind_speed_stacked, height_weights)
    A_z_ds = pd.DataFrame(
        {"A": A, "log_z": log_z}, index=wind_speed_stacked.x_y_time.to_index()
    ).to_xarray()

    return A_z_ds


def extrapolate(speeds, height_weights=None):
    """
    Get A and Z components of the log law describing wind speed
    at different heights above the ground, using recorded speed
    data at discrete heights above the ground.

    Uses sklearn LinearRegression function, find out more at:
    https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html

    Parameters
    ----------
    speeds : xr.DataArray
        dims = (height, all_other).
        If there are multiple non-height dimensions (e.g. lat, lon, time) then they stacked to a single dimension
        (i.e. `all_other` - name doesn't actually matter).
    height_weights : pd.Series or None (default).
        Weight per height. len(height_weights) = len(speeds.height)
        Weights force better fitting of curve to values of greater importance
        (e.g. to ensure that most likely wind turbine heights have closer matching values)
    """
    # It's a log law, so get some logging done
    log_height = np.log(speeds.height.values).reshape(-1, 1)

    if height_weights is not None:
        height_weights = height_weights.loc[speeds.height.values].values

    # Apply linear regression, which will apply the fit to every sample simultaneously (in parallel).
    # Weights give greater importance to fitting samples at certain heights.
    regression = LinearRegression().fit(
        log_height, speeds.values,
        sample_weight=height_weights
    )

    # extract coefficients
    # speed = A log(height) - A log(z),
    # so, A = coefficient (i.e. slope), z = exp(-intercept / A)
    A = regression.coef_.flatten()
    z = np.exp(np.divide(-regression.intercept_, A.flatten()))

    # remove small and large numbers
    # A[A < 0] = 0  # TODO: work out why this was commented out and if it should be added back
    z[z > 1e10] = 1e10
    z[z < 1e-10] = 1e-10

    return A.astype(np.float32), np.log(z).astype(np.float32)


if __name__ == '__main__':
    wind_speed_to_coeffs(
        path_to_input=snakemake.input.wind_data,
        height_weights=snakemake.params.height_weights,
        x_name=snakemake.params.x_name,
        y_name=snakemake.params.y_name,
        path_to_output=snakemake.output[0]
    )
