from string import ascii_lowercase

import xarray as xr
import matplotlib.pyplot as plt


CMAP = "magma"
TIMESERIES_COLORS = {
    "NEWA": "red",
    "COSMO-REA2": "blue",
    "Regridded COSMO-REA2": "purple",
    "Measured": "black",
}


plt.rcParams.update({
    "svg.fonttype": "none"
})

YEARS = slice("2011", "2013")  # period where there's overlap between all datasets


def plot_median_hourly_cfs(
    path_to_martigny_data, path_to_st_brais_data, path_to_output
):
    cfs = {
        "martigny": xr.open_dataset(path_to_martigny_data).capacityfactor,
        "st_brais": xr.open_dataset(path_to_st_brais_data).capacityfactor
    }

    fig = plt.figure(figsize=(15, 15))

    g = plt.GridSpec(
        nrows=5, ncols=2,
        figure=fig, height_ratios=[1, 10, 1, 10, 2], hspace=0.1
    )
    axes = {}

    handles, labels = plot_median_cfs(cfs, g, axes)

    axes["legend"] = plt.subplot(g[-1, :], frameon=False)
    axes["legend"].axis("off")
    axes["legend"].legend(
        handles=handles, labels=labels,
        bbox_to_anchor=(0.5, 0.5), loc="center", frameon=False, ncol=4
    )

    fig.savefig(path_to_output, pad_inches=0, bbox_inches="tight")


def median_cf(site_cf, months):
    return (
        site_cf
        .where(site_cf.time.dt.month.isin(months))
        .sel(time=YEARS)
        .groupby("time.hour").median()
        .to_series()
        .sum(level="hour")
    )


def plot_median_cfs(site_cfs, g, axes):
    i = 0
    letter_idx = 0
    for site, site_cf in site_cfs.items():
        axes[f'{site}_title'] = plt.subplot(g[i, :], frameon=False)
        axes[f'{site}_title'].axis('off')
        axes[f'{site}_title'].set_title(
            f"{ascii_lowercase[letter_idx]}. {site}",
            fontweight='bold', loc='left', y=0.5
        )
        i += 1
        j = 0
        for season, months in {"winter": [1, 2, 12], "summer": [6, 7, 8]}.items():
            if j == 0:
                axes[f'{site}_{season}'] = plt.subplot(g[i, j])
                shareax = axes[f'{site}_{season}']
                axes[f'{site}_{season}'].set_ylabel("Median capacity factor")
            else:
                axes[f'{site}_{season}'] = plt.subplot(g[i, j], sharey=shareax)
            axes[f'{site}_{season}'].set_xlabel("Hour of the day")
            axes[f'{site}_{season}'].legend().set_visible(False)
            for dataset in site_cf.dataset.values:
                median_ws = median_cf(site_cf.sel(dataset=dataset), months)
                median_ws.plot(
                    style=".",
                    color=TIMESERIES_COLORS[dataset],
                    label=dataset
                )
            if i == 1:
                axes[f'{site}_{season}'].set_title(season.title())
            handles, labels = axes[f'{site}_{season}'].get_legend_handles_labels()

            j += 1
        i += 1
        letter_idx += 1

    return handles, labels


if __name__ == "__main__":
    plot_median_hourly_cfs(
        path_to_martigny_data=snakemake.input.martigny_data,
        path_to_st_brais_data=snakemake.input.st_brais_data,
        path_to_output=snakemake.output[0]
    )
