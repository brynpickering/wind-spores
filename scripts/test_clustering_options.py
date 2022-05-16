from multiprocessing import Pool
from functools import partial

import pandas as pd

from cluster_cf_to_model_regions import (
    get_eligible_area_cf_ds,
    cluster_cf_data,
    get_region_cf_df
)


def plot_clusters(
    path_to_eligible_areas, path_to_turbine_cf,
    clustering_method, threads, path_to_output
):
    eligible_area_cf_ds = get_eligible_area_cf_ds(
        path_to_eligible_areas, path_to_turbine_cf
    )

    region_cf_df = get_region_cf_df(eligible_area_cf_ds, None)
    normalised_region_cf_df = region_cf_df.div(region_cf_df.mean())
    if clustering_method != "timeseries_kmeans":
        normalised_region_cf_df = normalised_region_cf_df.T

    cluster_numbers = list(range(5, 35, 5))

    with Pool(threads) as p:
        _func = partial(
            pooling_func,
            region_cf_df=normalised_region_cf_df,
            clustering_method=clustering_method
        )
        labels = p.map(_func, cluster_numbers)

    label_df = pd.DataFrame(
        index=cluster_numbers,
        data=labels,
        columns=region_cf_df.columns
    ).T
    label_df.to_csv(path_to_output)


def pooling_func(n_clusters, region_cf_df, clustering_method):
    return cluster_cf_data(region_cf_df, clustering_method, n_clusters)


if __name__ == "__main__":
    plot_clusters(
        path_to_turbine_cf=snakemake.input.turbine_cf,
        path_to_eligible_areas=snakemake.input.eligible_areas,
        clustering_method=snakemake.wildcards.clustering_method,
        threads=snakemake.threads,
        path_to_output=snakemake.output[0]
    )
