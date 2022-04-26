"""This module determines an upper bound of land eligibility for renewable generation based on geospatial data.

In here, we only exclude areas based on technical restrictions.
"""

import numpy as np
import rasterio

from renewablepotentialslib.eligibility import (
    FARM,
    FOREST,
    OTHER,
    Eligibility,
)

DATATYPE = np.float32


def determine_eligibility(
    path_to_land_cover, path_to_slope,
    path_to_building_share, path_to_urban_green_share, path_to_result,
    max_building_share, max_urban_green_share, tech
):
    """Determines eligibility of land for renewables."""
    with rasterio.open(path_to_land_cover) as src:
        transform = src.transform
        land_cover = src.read(1)
        crs = src.crs
    with rasterio.open(path_to_slope) as src:
        slope = src.read(1)
    with rasterio.open(path_to_building_share) as src:
        building_share = src.read(1)
    with rasterio.open(path_to_urban_green_share) as src:
        urban_green_share = src.read(1)
    eligibile_share = eligibility_land_mask(
        land_cover=land_cover,
        slope=slope,
        building_share=building_share,
        urban_green_share=urban_green_share,
        max_building_share=max_building_share,
        max_urban_green_share=max_urban_green_share,
        tech=tech
    )
    with rasterio.open(path_to_result, 'w', driver='GTiff', height=eligibile_share.shape[0],
                       width=eligibile_share.shape[1], count=1, dtype=DATATYPE,
                       crs=crs, transform=transform) as new_geotiff:
        new_geotiff.write(eligibile_share, 1)


def eligibility_land_mask(
    land_cover, slope, building_share, urban_green_share,
    max_building_share, max_urban_green_share, tech
):

    # prepare masks
    settlements = (building_share > max_building_share) | (urban_green_share > max_urban_green_share)
    farm = np.isin(land_cover, FARM)
    forest = np.isin(land_cover, FOREST)
    other = np.isin(land_cover, OTHER)
    if tech == "pv":
        slope[settlements | ~(farm | other)] = 0
    elif tech == "wind":
        slope[settlements | ~(farm | forest | other)] = 0
    return slope


if __name__ == "__main__":
    determine_eligibility(
        path_to_land_cover=snakemake.input.land_cover,
        path_to_slope=snakemake.input.slope,
        path_to_building_share=snakemake.input.building_share,
        path_to_urban_green_share=snakemake.input.urban_green_share,
        max_building_share=snakemake.params.max_building_share,
        max_urban_green_share=snakemake.params.max_urban_green_share,
        tech=snakemake.wildcards.tech,
        path_to_result=snakemake.output[0]
    )
