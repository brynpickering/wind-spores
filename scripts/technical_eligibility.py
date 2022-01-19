"""This module determines an upper bound of land eligibility for renewable generation based on geospatial data.

In here, we only exclude areas based on technical restrictions.
"""

import numpy as np
import rasterio

from renewablepotentialslib.eligibility import (
    DATATYPE,
    FARM,
    FOREST,
    OTHER,
    Eligibility,
    _add_eligibility
)


def determine_eligibility(
    path_to_land_cover, path_to_slope_pv, path_to_slope_wind,
    path_to_building_share, path_to_urban_green_share, path_to_result,
    max_building_share, max_urban_green_share, slope_threshold
):
    """Determines eligibility of land for renewables."""
    with rasterio.open(path_to_land_cover) as src:
        transform = src.transform
        land_cover = src.read(1)
        crs = src.crs
    with rasterio.open(path_to_slope_pv) as src:
        slope_pv = src.read(1)
    with rasterio.open(path_to_slope_wind) as src:
        slope_wind = src.read(1)
    with rasterio.open(path_to_building_share) as src:
        building_share = src.read(1)
    with rasterio.open(path_to_urban_green_share) as src:
        urban_green_share = src.read(1)
    eligibility = eligibility_land_mask(
        land_cover=land_cover,
        slope_pv=slope_pv,
        slope_wind=slope_wind,
        building_share=building_share,
        urban_green_share=urban_green_share,
        max_building_share=max_building_share,
        max_urban_green_share=max_urban_green_share,
        slope_threshold=slope_threshold
    )
    with rasterio.open(path_to_result, 'w', driver='GTiff', height=eligibility.shape[0],
                       width=eligibility.shape[1], count=1, dtype=DATATYPE,
                       crs=crs, transform=transform) as new_geotiff:
        new_geotiff.write(eligibility, 1)


def eligibility_land_mask(
    land_cover, slope_pv, slope_wind, building_share, urban_green_share,
    max_building_share, max_urban_green_share, slope_threshold
):
    # parameters
    assert (slope_pv <= slope_wind).all().all()  # wind can be built whereever pv can be built

    # prepare masks
    settlements = (building_share > max_building_share) | (urban_green_share > max_urban_green_share)
    farm = np.isin(land_cover, FARM)
    forest = np.isin(land_cover, FOREST)
    other = np.isin(land_cover, OTHER)
    pv = (slope_pv >= slope_threshold) & ~settlements & (farm | other)
    wind = (slope_wind >= slope_threshold) & ~settlements & (farm | forest | other)

    # allocate eligibility
    land = np.ones_like(land_cover, dtype=DATATYPE) * Eligibility.NOT_ELIGIBLE
    _add_eligibility(land, Eligibility.ONSHORE_WIND_AND_PV, wind & pv)
    _add_eligibility(land, Eligibility.ONSHORE_WIND, wind & ~pv)
    return land


if __name__ == "__main__":
    determine_eligibility(
        path_to_land_cover=snakemake.input.land_cover,
        path_to_slope_pv=snakemake.input.slope_pv,
        path_to_slope_wind=snakemake.input.slope_wind,
        path_to_building_share=snakemake.input.building_share,
        path_to_urban_green_share=snakemake.input.urban_green_share,
        max_building_share=snakemake.params.max_building_share,
        max_urban_green_share=snakemake.params.max_urban_green_share,
        slope_threshold=snakemake.params.slope_threshold,
        path_to_result=snakemake.output[0]
    )
