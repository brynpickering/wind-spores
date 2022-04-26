import xarray as xr
import pandas as pd


def regrid_dataset(
    path_to_timeseries_data, path_to_eligible_areas,
    dataset_config, level, cf_or_mwh, path_to_output
):
    data = xr.open_dataset(path_to_timeseries_data)
    area_df = pd.read_csv(path_to_eligible_areas)
    index_levels = [dataset_config["x-name"], dataset_config["y-name"]]
    ch_level = f"GID_{level}"

    weighted_ds = xr.concat(
        [regional_data(data, index_levels, area_df, ch_level, region, cf_or_mwh)
         for region in area_df[ch_level].unique()],
        dim="id"
    )
    weighted_ds.to_netcdf(path_to_output)


def regional_data(data, index_levels, area_df, ch_level, region, cf_or_mwh):
    area_weights = region_area_weights(index_levels, area_df, ch_level, region, cf_or_mwh)
    return (
        (data.interp(area_weights.coords) * area_weights)
        .sum(index_levels)
        .expand_dims({"id": [region]})
    )


def region_area_weights(index_levels, area_df, level, region, cf_or_mwh):
    # Get eligble area for all gridded data points in a particular region
    # return as percentage of total eligible area (CF) or total area (MWh)
    ds = (
        area_df[area_df[level] == region]
        .set_index(index_levels)
        .to_xarray()
    )
    if cf_or_mwh == "MWh":
        total_area_source = area_df[area_df[level] == region]
    elif cf_or_mwh == "CF":
        total_area_source = ds
    return ds.area / total_area_source.area.sum()


if __name__ == "__main__":
    regrid_dataset(
        path_to_timeseries_data=snakemake.input.timeseries_data,
        path_to_eligible_areas=snakemake.input.eligible_areas,
        dataset_config=snakemake.params.dataset_config,
        level=snakemake.wildcards.level,
        cf_or_mwh=snakemake.wildcards.cf_or_mwh,
        path_to_output=snakemake.output[0]
    )
