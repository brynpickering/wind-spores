import xarray as xr
import geopandas as gpd

import geoutil

CRS = "EPSG:3035"


def regrid_dataset(
    path_to_timeseries_data, path_to_polygons, path_to_eligible_land, path_to_ch_shape,
    dataset_config, level, cf_or_mwh, path_to_output
):
    data = xr.open_dataset(path_to_timeseries_data)

    polys_m = geoutil.load_polygons(path_to_polygons, dataset_config, CRS)

    ch_level = f"GID_{level}"
    ch_shape = (
        gpd.read_file(f"zip://{path_to_ch_shape}!gadm36_CHE.gpkg", layer=int(level))
        .to_crs(CRS)
        [[ch_level, "geometry"]]
    )
    polys_cut_to_ch_shape = gpd.overlay(polys_m.reset_index(), ch_shape)
    polys_cut_to_ch_shape = geoutil.get_eligible_land(polys_cut_to_ch_shape, path_to_eligible_land)

    polys_cut_to_ch_shape["area"] = polys_cut_to_ch_shape.area * polys_cut_to_ch_shape["fraction_eligible"]

    weighted_ds = xr.concat(
        [regional_data(data, polys_m, polys_cut_to_ch_shape, ch_level, region, cf_or_mwh)
         for region in ch_shape[ch_level].values],
        dim="id"
    )
    weighted_ds.to_netcdf(path_to_output)


def regional_data(data, polys, area_df, ch_level, region, cf_or_mwh):
    area_weights = region_area_weights(polys, area_df, ch_level, region, cf_or_mwh)
    return (
        (data.loc[area_weights.coords] * area_weights)
        .sum(polys.index.names)
        .expand_dims({"id": [region]})
    )


def region_area_weights(polys, area_df, level, region, cf_or_mwh):
    # Get eligble area for all gridded data points in a particular region
    # return as percentage of total eligible area (CF) or total area (MWh)
    ds = (
        area_df[area_df[level] == region]
        .set_index(polys.index.names)
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
        path_to_polygons=snakemake.input.polygons,
        path_to_eligible_land=snakemake.input.eligible_land,
        path_to_ch_shape=snakemake.input.ch_shape,
        dataset_config=snakemake.params.dataset_config,
        level=snakemake.wildcards.level,
        cf_or_mwh=snakemake.wildcards.cf_or_mwh,
        path_to_output=snakemake.output[0]
    )
