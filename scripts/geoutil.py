import geopandas as gpd
import numpy as np

from renewablepotentialslib.eligibility import Eligibility


def load_polygons(path_to_polygons, config, new_crs):
    polys = gpd.read_file(path_to_polygons, header=0)
    polys = polys.set_index([config["x-name"], config["y-name"]])
    polys = polys.to_crs(new_crs)
    return polys


def get_eligible_land(polygons, path_to_eligible_land):
    import rasterio
    import rasterio.mask

    eligible_wind_land = [Eligibility.ONSHORE_WIND, Eligibility.ONSHORE_WIND_AND_PV]
    with rasterio.open(path_to_eligible_land, "r") as src:
        _meta = src.meta
        polygons_to_src_crs = polygons.to_crs(_meta["crs"])
        for idx in polygons.index:
            geom = polygons_to_src_crs.loc[idx, "geometry"]
            _out, _ = rasterio.mask.mask(src, [geom], crop=True, nodata=0)
            fraction_eligible = np.isin(_out, eligible_wind_land).sum() / _out.size
            polygons.loc[idx, "fraction_eligible"] = fraction_eligible
    return polygons
