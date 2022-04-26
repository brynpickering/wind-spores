import pandas as pd
import seaborn as sns
import numpy as np

import xarray as xr
import matplotlib.pyplot as plt

COLOURS = {
    "Enercon E53 800": "#fde725",
    "enercon_e53_800_60": "#fde725",
    "Enercon E82 3000": "#90d743",
    "enercon_e82_3000_60": "#90d743",
    "enercon_e82_3000_80": "#35b779",
    "Vestas V112 3450": "#21918c",
    "vestas_v112_3450_80": "#21918c",
    "Vestas V110 2000": "#31688e",
    "vestas_v110_2000_80": "#31688e",
    "vestas_v110_2000_100": "#443983",
    "vestas_v110_2000_120": "#440154",
}

plt.rcParams.update({
    "svg.fonttype": "none"
})


def plot_study_turbine_comparison(paths_to_cfs, turbine_config, path_to_power_curves, path_to_output):
    capacityfactor_dataset = xr.concat(
        [xr.open_dataset(file) for file in paths_to_cfs], dim="turbine"
    )

    capacityfactor_df = (
        capacityfactor_dataset
        .capacityfactor
        .to_series()
        .unstack("turbine")
        .droplevel("id")
        .reindex(COLOURS.keys(), axis=1)
        .dropna(how="all", axis=1)
    )
    legend_labels = [
        f"{turbine_config[turbine]['long-name']} @ {height}m"
        for turbine, height
        in capacityfactor_df.columns.str.rsplit("_", 1)
    ]
    power_curves = pd.read_csv(path_to_power_curves, index_col=0)

    with sns.plotting_context("paper", font_scale=1.18):
        fig = plot_all_panels(capacityfactor_df, legend_labels, power_curves)

    fig.savefig(path_to_output, pad_inches=0, bbox_inches="tight")


def plot_all_panels(capacityfactor_df, legend_labels, power_curves):

    fig = plt.figure(figsize=(18, 15))
    g = plt.GridSpec(
        nrows=2, ncols=2,
        figure=fig, hspace=0.15
    )
    axes = {}

    axes["power_curves"] = plt.subplot(g[0, 0], frameon=False)
    axes["load_duration"] = plt.subplot(g[0, 1], frameon=False)
    axes["cf_boxplots"] = plt.subplot(g[1, :], frameon=False)

    plot_load_duration_curve(capacityfactor_df, axes["load_duration"], legend_labels)
    plot_cf_boxplots(capacityfactor_df, axes["cf_boxplots"], legend_labels)
    plot_power_curves(power_curves, axes["power_curves"])

    titles = {
        "power_curves": "a. Turbine power curves",
        "load_duration": "b. Swiss average capacity factor duration curve",
        "cf_boxplots": "c. Swiss average hourly capacity factor distribution",
    }
    for ax_name, title in titles.items():
        axes[ax_name].set_title(title, fontweight='bold', loc='left', y=1)

    return fig


def plot_load_duration_curve(capacityfactor_df, ax, legend_labels):
    df = pd.concat(
        [capacityfactor_df[i].sort_values(ascending=False).reset_index(drop=True)
         for i in capacityfactor_df.columns],
        keys=capacityfactor_df.columns, axis=1
    )
    df.plot(color=[COLOURS[i] for i in df.columns], lw=2, ax=ax)
    sns.despine(ax=ax)
    ax.set_xlabel("Hour of year")
    ax.set_ylabel("Capacity factor")
    ax.legend(frameon=False)
    ax.set_ylim((0, 1.05))
    ax.set_xticks(np.linspace(0, 61320, 12))
    ax.set_xticklabels([0] + ["" for i in range(10)] + [8760])
    ax.legend(labels=legend_labels, frameon=False)


def _get_season(x):
    if x.month in [1, 2, 12]:
        return "winter"
    elif x.month in [6, 7, 8]:
        return "summer"
    else:
        return "shoulder"


def plot_cf_boxplots(capacityfactor_df, ax, legend_labels):
    df_violin = (
        capacityfactor_df
        .groupby([capacityfactor_df.index, _get_season]).mean()
        .drop("shoulder", level=1)
        .rename_axis(index=["time", "season"])
        .stack()
        .to_frame("CF")
        .reset_index()
    )
    df_box = (
        capacityfactor_df
        .stack()
        .to_frame("CF")
        .reset_index()
    )
    df_median = (
        df_box
        .groupby("turbine")
        .median()
        .reindex(capacityfactor_df.columns)
        .reset_index()
    )

    sns.violinplot(
        data=df_violin,
        x="turbine",
        y="CF",
        hue="season",
        split=True,
        scale="area",
        cut=0,
        ax=ax,
        inner=None,
        zorder=-1,
        order=capacityfactor_df.columns
    )
    sns.boxplot(
        data=df_box,
        x="turbine",
        y="CF",
        zorder=1,
        width=0.08,
        palette=COLOURS,
        fliersize=0,
        linewidth=0,
        order=capacityfactor_df.columns,
        ax=ax,
        boxprops={'zorder': 1, "lw": 0.3, "ec": "white"},
        whiskerprops={'zorder': 1, "color": "white", "lw": 1},
        capprops={'zorder': 1, "color": "white", "lw": 1},
        medianprops={'lw': 0, "color": "white"}
    )
    sns.scatterplot(
        data=df_median,
        x="turbine",
        y="CF",
        ax=ax,
        color="white",
        lw=0,
        s=5,
        zorder=2
    )

    sns.despine(ax=ax)
    ax.set_xlabel("Turbine")
    ax.set_xticklabels(legend_labels)
    ax.set_ylabel("Capacity factor")
    ax.legend(frameon=False, ncol=2, title="Season", bbox_to_anchor=(0.5, 1), loc="center")
    ax.set_ylim((0, 1))


def plot_power_curves(power_curves, ax):
    df = (
        power_curves
        .reindex(COLOURS.keys(), axis=1)
        .dropna(how="all", axis=1)
        .reindex(np.linspace(0, 35, 351))
        .interpolate()
    )
    df = df.where((df > 0) & (df.diff() >= 0))
    df.plot(
        color=[COLOURS[i] for i in df.columns], lw=2, ax=ax
    )
    sns.despine(ax=ax)
    ax.set_xlabel("Wind speed (m/s)")
    ax.set_ylabel("Capacity factor")
    ax.legend(frameon=False)
    ax.set_ylim((0, 1.05))


if __name__ == "__main__":
    plot_study_turbine_comparison(
        paths_to_cfs=snakemake.input.cfs,
        path_to_power_curves=snakemake.input.power_curves,
        turbine_config=snakemake.params.turbine_config,
        path_to_output=snakemake.output[0]
    )
