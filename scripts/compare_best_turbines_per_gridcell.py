import logging

import pandas as pd
import xarray as xr
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import geoutil

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S'
)


CMAP = "magma"
CMAP_CORR = "PuOr"
TIMESERIES_COLORS = {
    "newa": "red",
    "cosmo-rea2": "blue",
}
TITLES = {
    'q0.1': "Max 10th percentile",
    'q0.25': "Max 25th percentile",
    'mean': "Max mean",
    'q0.75': "Max 75th percentile",
    'q0.9': "Max 90th percentile",
    'q1.0': "Max max",
    'std': "Min standard deviation"
}


plt.rcParams.update({
    "svg.fonttype": "none"
})


def compare_turbines(
    path_to_top_turbines, path_to_ch_shape, path_to_polys,
    plot_crs, dataset_config, turbine_config, level, path_to_output
):
    top_turbines = xr.open_dataset(path_to_top_turbines)
    ch_shape = gpd.read_file(f"zip://{path_to_ch_shape}!gadm36_CHE.gpkg", layer=0).to_crs(plot_crs)
    if level is None:
        polys = geoutil.load_polygons(path_to_polys, dataset_config, plot_crs)
    elif path_to_polys is None:
        polys = (
            gpd.read_file(f"zip://{path_to_ch_shape}!gadm36_CHE.gpkg", layer=int(level))
            .to_crs(plot_crs)
            .set_index(f"GID_{level}")
        )
    polys["ch"] = polys.intersects(ch_shape.geometry[0])

    fig = plt.figure(figsize=(15, 24))
    g = plt.GridSpec(
        9, 4,
        height_ratios=[1, 10, 10, 10, 10, 10, 10, 10, 4], width_ratios=[1, 10, 10, 10],
        hspace=0.05, wspace=0.05
    )
    ax = {}
    set_metric_titles(ax, g)
    plot_compare_turbines(ax, g, polys, turbine_config, ch_shape, top_turbines)

    fig.savefig(path_to_output, bbox_inches="tight")


def set_metric_titles(ax, g):
    _row = 1
    _column = 0
    for metric in ['q0.1', 'q0.25', 'mean', 'q0.75', 'q0.9', 'q1.0', 'std']:
        _ax = ax[f"{metric}_title"] = plt.subplot(g[_row, _column], frameon=False)
        _ax.annotate(TITLES[metric], xy=(0.5, 0.5), xycoords="axes fraction", rotation=90, ha="center", va="center", fontsize=16)
        _ax.axis("off")
        _row += 1


def plot_compare_turbines(ax, g, polys, turbine_config, ch_shape, top_turbines):

    def _turbine_name(short_name):
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
    colours = get_cmap(turbine_config, top_turbines)

    _row = 0
    _column = 1
    for season in ["all", "summer", "winter"]:
        _ax = ax[f"{season}_title"] = plt.subplot(g[_row, _column], frameon=False)
        _ax.annotate(season, xy=(0.5, 0.5), xycoords="axes fraction", ha="center", va="center", fontsize=16)
        _ax.axis("off")
        _row = 1
        for metric in ['q0.1', 'q0.25', 'mean', 'q0.75', 'q0.9', 'q1.0', 'std']:
            _ax = ax[f"{season}_{metric}"] = plt.subplot(g[_row, _column], frameon=False)
            polys = add_best_to_polys(polys, metric, season, top_turbines)
            valid_polys = polys[(polys.best_turbine != "")]
            _c = []
            for i in valid_polys.best_turbine:
                _c.append(colours[i])
            valid_polys.plot(
                color=_c,
                legend=False,
                antialiased=False,
                ax=_ax
            )
            _ax.annotate(
                _turbine_name(
                    valid_polys[valid_polys.ch]
                    .best_turbine
                    .value_counts()
                    .index[0]
                ),
                xy=(0, 0.95), xycoords="axes fraction", fontsize=10,
                bbox={"facecolor": "white", "alpha": 0.5, "pad": 0, "edgecolor": "None"}
            )
            xlim = [ch_shape.total_bounds[0] * 0.999, ch_shape.total_bounds[2] * 1.001]
            ylim = [ch_shape.total_bounds[1] * 0.999, ch_shape.total_bounds[3] * 1.001]
            _ax.set_xlim(*xlim)
            _ax.set_ylim(*ylim)
            ch_shape.plot(ax=_ax, ec="black", fc="None")
            _ax.axis("off")
            if _row == 7:
                _row = 0
                _column += 1
            else:
                _row += 1
    add_legend(ax, g, colours, _turbine_name)


def add_legend(ax, g, colours, turbine_namer):
    _ax = ax["legend"] = plt.subplot(g[-1, :], frameon=False)
    _ax.axis("off")
    patches = [
        mpl.patches.Patch(facecolor=_color, label=turbine_namer(_turbine))
        for _turbine, _color in colours.items()
    ]
    _ax.legend(
        handles=patches,
        ncol=int(len(colours) / 3),
        bbox_to_anchor=(0.5, 0.5),
        loc="center",
        frameon=False
    )


def get_cmap(turbine_config, top_turbines):
    colours = {}

    # study turbine config doesn't differentiate by turbine height
    # so we link the 'top turbines' which have heights in their names
    # with the relevant config item
    unique_config_items = [
        i for i in np.unique(top_turbines.to_array()) if i != ""
    ]
    if len(unique_config_items) > len(turbine_config.keys()):
        colour_keys = []
        for i in turbine_config.keys():
            colour_keys.extend([
                j for j in sorted(unique_config_items)
                if i in j
            ])
    else:
        colour_keys = turbine_config.keys()

    palette = sns.color_palette("magma", n_colors=len(colour_keys))
    i = 0
    for k in colour_keys:
        colours[k] = palette[i]
        i += 1
    return colours


def add_best_to_polys(polys, metric, period, top_turbines):
    if metric == "std":
        var = "_min"
    else:
        var = "_max"
    best_turbine = (
        top_turbines
        .sel(metric=metric, period=period)
        [f"turbine{var}"]
        .to_series()
    )
    if isinstance(best_turbine.index, pd.MultiIndex):
        best_turbine = best_turbine.reorder_levels(polys.index.names)

    polys["best_turbine"] = best_turbine.reindex(polys.index)
    return polys


if __name__ == "__main__":
    compare_turbines(
        path_to_top_turbines=snakemake.input.top_turbines,
        path_to_ch_shape=snakemake.input.ch_shape,
        path_to_polys=snakemake.input.get("polys", None),
        level=snakemake.params.level,
        plot_crs=snakemake.params.plot_crs,
        dataset_config=snakemake.params.dataset_config,
        turbine_config=snakemake.params.turbine_config,
        path_to_output=snakemake.output[0]
    )
