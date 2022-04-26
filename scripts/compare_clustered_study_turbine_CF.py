import pandas as pd
import geopandas as gpd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import geoutil


def plot_turbine_clustered_cf_comparison(
    path_to_polygons, path_to_eligible_areas, path_to_model_region_mapping, path_to_ch_shape,
    paths_to_cf, paths_to_cluster_labels, plot_crs, turbine_config, dataset_config, path_to_output
):

    polygons = geoutil.load_polygons(path_to_polygons, dataset_config, plot_crs)
    polygon_da = polygons.to_xarray()

    available_area_da = (
        pd.read_csv(path_to_eligible_areas, index_col=["rlat", "rlon", "GID_1"])
        .to_xarray()
        .interp(rlat=polygon_da.rlat, rlon=polygon_da.rlon, method="nearest")
    )
    available_area_da["geometry"] = polygon_da.geometry
    cf_timeseries = {
        get_turbine_name(path_to_cf): pd.read_csv(
            path_to_cf, index_col=0, parse_dates=True, header=[0, 1]
        )
        for path_to_cf in paths_to_cf
    }
    region_mapping = pd.read_csv(path_to_model_region_mapping)
    ch_canton_shapes = geoutil.load_ch_shapes(path_to_ch_shape, level=1).to_crs(plot_crs)
    region_shapes = (
        ch_canton_shapes
        .replace(region_mapping.set_index("gadm").model_region.to_dict())
        .dissolve("GID_1")
    )

    with sns.plotting_context("paper", font_scale=1.2):
        ncols = len(paths_to_cf)
        fig, ax = plt.subplots(
            4, ncols, figsize=(4 * ncols, 10),
            gridspec_kw={"height_ratios": [1, 20, 1, 30], "hspace": 0.1, "wspace": 0.1}
        )

        ax[0, 0].annotate(
            "a. average capacity factor per cluster",
            xy=(0, 0), xycoords="axes fraction", ha="left", va="bottom", fontsize=12,
            fontweight="bold"
        )
        ax[0, 0].axis("off")
        _col = 0
        for path_to_cluster_labels in paths_to_cluster_labels:
            ax[0, _col].axis("off")
            turbine_name = get_turbine_name(path_to_cluster_labels)
            plot_maps(
                ax[1, _col], cf_timeseries,
                available_area_da, polygon_da, region_shapes, path_to_cluster_labels,
                turbine_name, turbine_config, plot_crs
            )
            _col += 1

        ax[2, 0].annotate(
            "b. hourly capacity factor duration curves in every cluster",
            xy=(0, -0.2), xycoords="axes fraction", ha="left", va="bottom", fontsize=12,
            fontweight="bold"
        )

        _col = 0
        for turbine_name, turbine_cf in cf_timeseries.items():
            ax[2, _col].axis("off")
            if _col == 0:
                label_yaxis = True
            else:
                label_yaxis = False
            plot_duration_curves(ax[3, _col], turbine_cf, label_yaxis)
            _col += 1

        fig.savefig(path_to_output, bbox_inches="tight", dpi=300)


def get_turbine_name(path_to_file):
    return path_to_file.split("/")[-1].split("-")[0]


def plot_maps(
    ax, cf_dfs, available_area_da, polygon_da, region_shapes, path_to_cluster_labels,
    turbine, turbine_config, plot_crs
):
    available_area_da["clusters"] = (
        pd.read_csv(path_to_cluster_labels)
        .set_index(["rlat", "rlon", "id"])
        .squeeze()
        .rename_axis(index={"id": "GID_1"})
        .to_xarray()
        .interp(rlat=polygon_da.rlat, rlon=polygon_da.rlon, method="nearest")
    )
    cluster_gdf = gpd.GeoDataFrame(
        available_area_da.to_dataframe().dropna(subset=["area"]).reset_index(),
        crs=plot_crs
    )

    turbine_cluster_gdf = cluster_gdf.dissolve(["GID_1", "clusters"])
    mean_cf = cf_dfs[turbine].mean().rename_axis(index={"id": "GID_1", "cluster": "clusters"})
    mean_cf.index = mean_cf.index.set_levels(mean_cf.index.levels[1].astype(float), level=1)
    turbine_cluster_gdf["mean"] = mean_cf.dropna().reindex(turbine_cluster_gdf.index)

    turbine_cluster_gdf.plot(
        "mean", cmap="viridis", antialiased=False, ax=ax, legend=True,
        legend_kwds={'orientation': 'horizontal', 'label': "Capacity factor", 'shrink': 0.8}
    )
    region_shapes.plot(fc="none", ec="red", lw=1, linestyle="--", ax=ax)
    ax.axis("off")
    ax.set_title(get_turbine_long_name(turbine, turbine_config))


def plot_duration_curves(ax, cf_df, label_yaxis):
    ordered_cf_df = pd.concat([
        cf_df.loc[:, i].sort_values(ascending=False).reset_index(drop=True)
        for i in cf_df.columns
    ], axis=1)
    ordered_cf_df.plot(color="blue", alpha=0.3, legend=False, ax=ax)
    if label_yaxis:
        ax.set_ylabel("Capacity factor")
    else:
        ax.set_ylabel("")
        ax.set_yticklabels([])
    ax.set_xlabel("Hour of year")
    sns.despine(ax=ax)
    ax.set_ylim((0, 1.05))
    ax.set_xticks(np.linspace(0, 61320, 12))
    ax.set_xticklabels([0] + ["" for i in range(10)] + [8760])


def get_turbine_long_name(short_name, turbine_config):
    if short_name not in turbine_config.keys():
        short_name, height = short_name.rsplit("_", 1)
        return (
            f"{turbine_config[short_name]['long-name']} @ "
            f"{height}m"
        )
    else:
        return (
            f"{turbine_config[short_name]['long-name']} @ "
            f"{turbine_config[short_name]['height']}m"
        )

if __name__ == "__main__":
    plot_turbine_clustered_cf_comparison(
        path_to_polygons=snakemake.input.polygons,
        path_to_eligible_areas=snakemake.input.eligible_areas,
        path_to_model_region_mapping=snakemake.input.model_region_mapping,
        path_to_ch_shape=snakemake.input.ch_shape,
        paths_to_cf=snakemake.input.cf,
        paths_to_cluster_labels=snakemake.input.cluster_labels,
        plot_crs=snakemake.params.plot_crs,
        turbine_config=snakemake.params.turbine_config,
        dataset_config=snakemake.params.dataset_config,
        path_to_output=snakemake.output[0]
    )
