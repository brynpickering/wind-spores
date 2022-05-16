import geopandas as gpd


def load_polygons(path_to_polygons, config, new_crs):
    polys = gpd.read_file(path_to_polygons, header=0)
    polys = polys.set_index([config["x-name"], config["y-name"]])
    polys = polys.to_crs(new_crs)
    return polys


def get_eligible_land(polygons, path_to_eligible_land):
    import rasterio
    import rasterio.mask
    from rasterstats import zonal_stats

    with rasterio.open(path_to_eligible_land, "r") as src:
        _meta = src.meta
        polygons_to_src_crs = polygons.to_crs(_meta["crs"])
        potentials = zonal_stats(
            polygons_to_src_crs,
            src.read(1),
            affine=_meta["transform"],
            stats="mean",
            nodata=-999
        )
        average_land_area_share = [stat["mean"] for stat in potentials]
        polygons["fraction_eligible"] = average_land_area_share

    return polygons


def load_ch_shapes(path_to_ch_shapes, level):
    return gpd.read_file(f"zip://{path_to_ch_shapes}!gadm36_CHE.gpkg", layer=int(level))
