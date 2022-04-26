import geopandas as gpd
import pandas as pd

import geoutil

CRS = "EPSG:3035"


def get_eligible_area_per_gricell(
    path_to_polygons, path_to_ch_shape, path_to_eligible_land,
    level, dataset_config, path_to_output
):
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

    polys_cut_to_ch_shape.set_index(["rlon", "rlat", "GID_1"])["area"].to_csv(path_to_output)


def ch_level_eligible_area_to_model_regions(
    path_to_model_region_mapping, path_to_eligible_areas, path_to_output
):
    region_mapping_df = pd.read_csv(path_to_model_region_mapping)

    area_df = pd.read_csv(path_to_eligible_areas)

    # map from CHE cantons to model regions
    renamed_area_df = area_df.replace(region_mapping_df.set_index("gadm").model_region.to_dict())

    model_area_df = (
        renamed_area_df
        .where(renamed_area_df.area > 0)
        .dropna()
        .groupby(["rlat", "rlon", "GID_1"]).sum()
    )
    model_area_df.to_csv(path_to_output)


if __name__ == "__main__":
    if "model_region_mapping" in snakemake.input.keys():
        ch_level_eligible_area_to_model_regions(
            path_to_eligible_areas=snakemake.input.eligible_areas,
            path_to_model_region_mapping=snakemake.input.model_region_mapping,
            path_to_output=snakemake.output[0]
        )
    else:
        get_eligible_area_per_gricell(
            path_to_polygons=snakemake.input.polygons,
            path_to_ch_shape=snakemake.input.ch_shape,
            path_to_eligible_land=snakemake.input.eligible_land,
            level=snakemake.wildcards.level,
            dataset_config=snakemake.params.dataset_config,
            path_to_output=snakemake.output[0]
        )

