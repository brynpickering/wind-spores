import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt

import geoutil


def plot_clusters(
    path_to_model_region_mapping, path_to_label_df,
    path_to_polygons, path_to_ch_shape,
    dataset_config, plot_crs, clustering_method,
    path_to_output,
):
    idx_items = [dataset_config["x-name"], dataset_config["y-name"]]
    label_df = pd.read_csv(path_to_label_df, index_col=idx_items, header=0)
    label_df.columns = label_df.columns.astype(int)

    model_region_shapes, label_shapes = get_geometries(
        path_to_polygons, dataset_config,
        path_to_ch_shape, path_to_model_region_mapping,
        label_df, plot_crs
    )

    plot_geometries(label_shapes, model_region_shapes, clustering_method, path_to_output)


def get_geometries(
    path_to_polygons, dataset_config, path_to_ch_shape, path_to_model_region_mapping,
    label_df, plot_crs
):
    polys_m = geoutil.load_polygons(path_to_polygons, dataset_config, plot_crs)
    ch_level = "GID_1"
    ch_shape = (
        gpd.read_file(f"zip://{path_to_ch_shape}!gadm36_CHE.gpkg", layer=1)
        .to_crs(plot_crs)
        [[ch_level, "geometry"]]
    )
    region_mapping_df = pd.read_csv(path_to_model_region_mapping)
    ch_shape["model_region"] = (
        region_mapping_df
        .set_index("gadm")
        .reindex(ch_shape.GID_1.values)
        .model_region
        .values
    )

    model_region_shapes = ch_shape.dissolve("model_region")

    label_ds = (
        label_df
        .to_xarray()
        .interp(rlat=polys_m.to_xarray().rlat, rlon=polys_m.to_xarray().rlon, method="nearest")
    )

    label_ds["geometry"] = polys_m.to_xarray().geometry

    label_shapes = gpd.GeoDataFrame(label_ds.to_dataframe(), crs=plot_crs)

    return model_region_shapes, label_shapes


def plot_geometries(label_shapes, model_region_shapes, clustering_method, path_to_output):
    fig, ax = plt.subplots(3, 2, figsize=(15, 15))
    fig.suptitle(f"Clustering using `{clustering_method}` algorithm")
    _row = 0
    _col = 0
    for i in [5, 10, 15, 20, 25, 30]:
        label_shapes.plot(i, ax=ax[_row, _col], antialiased=False)
        model_region_shapes.plot(ax=ax[_row, _col], fc="None", ec="red", linestyle="--")
        ax[_row, _col].set_title(f"{i} clusters")
        ax[_row, _col].axis("off")
        if _col == 1:
            _row += 1
            _col = 0
        else: _col += 1

    fig.savefig(path_to_output, bbox_inches="tight", dpi=300)


if __name__ == "__main__":
    plot_clusters(
        path_to_model_region_mapping=snakemake.input.model_region_mapping,
        path_to_label_df=snakemake.input.label_df,
        path_to_polygons=snakemake.input.polygons,
        path_to_ch_shape=snakemake.input.ch_shape,
        clustering_method=snakemake.wildcards.clustering_method,
        dataset_config=snakemake.params.dataset_config,
        plot_crs=snakemake.params.plot_crs,
        path_to_output=snakemake.output[0],
    )
