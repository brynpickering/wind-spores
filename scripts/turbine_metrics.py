import logging

import xarray as xr

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S'
)


def turbine_metric(
    path_to_turbine_dataset, comparison_quantiles, comparison_periods, path_to_output
):

    turbine_cf = xr.open_dataset(path_to_turbine_dataset).capacityfactor

    results = []
    for period, months in comparison_periods.items():
        logging.info(period)
        period_specific_turbine_da = (
            turbine_cf.where(turbine_cf.time.dt.month.isin(months))
        )
        logging.info("getting quantiles")
        period_specific_turbine_quantiles = (
            period_specific_turbine_da.quantile(comparison_quantiles, dim="time")
        )
        logging.info("getting mean")
        period_specific_turbine_mean = (
            period_specific_turbine_da.mean(dim="time").expand_dims({"quantile": ["mean"]})
        )

        logging.info("combining metrics")
        period_specific_turbine_metrics = (
            xr.concat(
                [period_specific_turbine_mean, period_specific_turbine_quantiles],
                dim="quantile"
            )
            .expand_dims({"period": [period]})
        )
        period_specific_turbine_metrics.coords["quantile"] = (
            period_specific_turbine_metrics.coords["quantile"].astype(str)
        )
        results.append(period_specific_turbine_metrics)
    logging.info("combining periods")
    result_da = xr.concat(results, dim="period")

    result_da.rename("turbine_metrics").to_netcdf(path_to_output)


if __name__ == "__main__":
    turbine_metric(
        path_to_turbine_dataset=snakemake.input.turbine_output,
        comparison_quantiles=snakemake.params.comparison_quantiles,
        comparison_periods=snakemake.params.comparison_periods,
        path_to_output=snakemake.output[0]
    )
