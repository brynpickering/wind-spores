import logging

import xarray as xr
import geopandas as gpd
import geoutil

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S'
)


def turbine_metric(
    path_to_turbine_dataset, path_to_polygons, path_to_eligible_land,
    mask_on_eligibility, comparison_quantiles, comparison_periods, turbine_capacity,
    cf_or_mwh, turbine_density, dataset_config, level, path_to_output
):
    if level is None:
        polys_m = geoutil.load_polygons(path_to_polygons, dataset_config, "EPSG:3035")
    else:
        polys_m = (
            gpd.read_file(f"zip://{path_to_polygons}!gadm36_CHE.gpkg", layer=int(level))
            .to_crs("EPSG:3035")
            .set_index(f"GID_{level}")
            .rename_axis(index="id")
        )
    turbine_output_per_gridcell = get_turbine_output_per_gridcell(
        path_to_turbine_dataset, turbine_capacity, polys_m, cf_or_mwh,
        mask_on_eligibility, path_to_eligible_land, dataset_config, turbine_density
    )

    results = []
    for period, months in comparison_periods.items():
        logging.info(period)
        period_specific_turbine_da = (
            turbine_output_per_gridcell
            .where(turbine_output_per_gridcell.time.dt.month.isin(months))
        )

        period_specific_turbine_metrics = (
            xr.concat(
                [
                    get_quantiles(period_specific_turbine_da, comparison_quantiles),
                    get_metric(period_specific_turbine_da, "mean"),
                    get_metric(period_specific_turbine_da, "std") / get_metric(period_specific_turbine_da, "mean", add_dim=False)  # normalised standard deviation
                ],
                dim="metric"
            )
            .expand_dims({"period": [period]})
        )
        results.append(period_specific_turbine_metrics)

    logging.info("combining periods")

    result_da = xr.concat(results, dim="period")

    # ignore gridcells without useful data
    result_da = result_da.where(turbine_output_per_gridcell.sum("time", min_count=1).notnull())

    result_da = result_da.assign_attrs(unit=cf_or_mwh, masked=mask_on_eligibility)
    result_da.rename("turbine_metrics").to_netcdf(path_to_output)


def get_quantiles(da, quantiles):
    logging.info("getting quantiles")
    da_quantiles = da.quantile(quantiles, dim="time")

    da_quantiles.coords["quantile"] = da_quantiles["quantile"].astype(str)

    da_quantiles.coords["quantile"] = "q" + da_quantiles["quantile"].to_index()

    return da_quantiles.rename({"quantile": "metric"})


def get_metric(da, metric, add_dim=True):
    logging.info(f"getting {metric}")
    da_metric = getattr(da, metric)(dim="time")
    if add_dim:
        return da_metric.expand_dims({"metric": [metric]})
    else:
        return da_metric


def get_turbine_output_per_gridcell(
    path_to_turbine_dataset, turbine_capacity, polys_m, cf_or_mwh,
    mask_on_eligibility, path_to_eligible_land, dataset_config, turbine_density
):
    turbine_cf = xr.open_dataset(path_to_turbine_dataset).capacityfactor
    if cf_or_mwh == "CF":
        return turbine_cf

    turbine_output = turbine_cf * turbine_capacity
    if mask_on_eligibility == "masked":
        polys_m["area"] = polys_m.area * geoutil.get_eligible_land(polys_m, path_to_eligible_land)["fraction_eligible"]
    else:
        polys_m["area"] = polys_m.area

    max_turbine_MW_per_gridcell = polys_m["area"].div(1000**2).mul(turbine_density).to_xarray()
    turbines_per_gridcell = max_turbine_MW_per_gridcell / turbine_capacity

    # We want to explicitly ignore gridcells with no turbines in downstream metrics
    turbines_per_gridcell = turbines_per_gridcell.where(turbines_per_gridcell > 0)

    return turbines_per_gridcell * turbine_output


if __name__ == "__main__":
    turbine_metric(
        path_to_turbine_dataset=snakemake.input.turbine_output,
        path_to_polygons=snakemake.input.polygons,
        path_to_eligible_land=snakemake.input.eligible_land,
        turbine_capacity=snakemake.params.turbine_capacity,
        turbine_density=snakemake.params.turbine_density,
        dataset_config=snakemake.params.dataset_config,
        comparison_quantiles=snakemake.params.comparison_quantiles,
        comparison_periods=snakemake.params.comparison_periods,
        level=snakemake.params.level,
        mask_on_eligibility=snakemake.wildcards.ismasked,
        cf_or_mwh=snakemake.wildcards.cf_or_mwh,
        path_to_output=snakemake.output[0]
    )
