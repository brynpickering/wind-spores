from multiprocessing import Pool
from functools import partial

import xarray as xr
import numpy as np

from sklearn.cluster import SpectralClustering, KMeans
#from tslearn.clustering import TimeSeriesKMeans

import pandas as pd

CLUSTERING_KWARGS = {
    "kmeans": {},
    "spectral": {"affinity": "nearest_neighbors", "assign_labels": "kmeans"},
    "timeseries_kmeans": {"metric": "dtw"}
}

def create_clusters(
    path_to_turbine_cf, path_to_eligible_areas,
    clustering_method, model_region_config, threads,
    path_to_timeseries_output, path_to_label_output
):


    eligible_area_cf_ds = get_eligible_area_cf_ds(
        path_to_eligible_areas, path_to_turbine_cf
    )

    with Pool(threads) as p:
        _func = partial(
            pooling_func,
            eligible_area_cf_ds=eligible_area_cf_ds,
            clustering_method=clustering_method
        )
        cluster_data = p.map(_func, list(model_region_config.items()))
    cluster_timeseries_df = pd.concat([i[0] for i in cluster_data], axis=1)
    label_df = pd.concat([i[1] for i in cluster_data])

    cluster_timeseries_df.to_csv(path_to_timeseries_output)
    label_df.to_csv(path_to_label_output)


def pooling_func(model_region_config, eligible_area_cf_ds, clustering_method):
    model_region, region_config = model_region_config
    n_clusters = region_config["max-clusters"]
    region_cf_df = get_region_cf_df(eligible_area_cf_ds, model_region)
    normalised_region_cf_df = region_cf_df.div(region_cf_df.mean())
    if clustering_method != "timeseries_kmeans":
        normalised_region_cf_df = normalised_region_cf_df.T
    labels = cluster_cf_data(normalised_region_cf_df, clustering_method, n_clusters)

    label_df = pd.Series(index=region_cf_df.columns, data=labels)
    label_df_to_return = (
        label_df
        .to_frame(model_region)
        .rename_axis(columns="id")
        .stack()
    )

    areas = get_region_area(eligible_area_cf_ds, model_region)
    return get_cluster_timeseries(label_df, region_cf_df, areas, model_region), label_df_to_return


def get_cluster_timeseries(label_df, region_cf_df, areas, region_name):
    weighted_average_cf = (
        region_cf_df
        .groupby(label_df, axis=1)
        .apply(get_weighted_average, weights=areas)
        .rename_axis(index="timesteps", columns="cluster")
    )
    weighted_average_cf = pd.concat(
        [weighted_average_cf], axis=1, keys=[region_name], names=["id"]
    )
    return weighted_average_cf


def get_eligible_area_cf_ds(path_to_eligible_areas, path_to_turbine_cf):
    """
    Align per-gridcell capacity factor data with per-gridcell available area for wind power
    deployment. This will remove gridcells outside switzerland, those entirely over lakes,
    and those in urban centres.

    Args:
        path_to_eligible_areas (str): path to CSV file containing a column for each
            gridcell x and y coordinate (e.g. "rlat" and "rlon"),
            "GID_1" (Euro-Calliope regions based on aggregated swiss GADM canton codes),
            and "area" (eligible area data).
        path_to_turbine_cf (str): path to xarray dataset with datarray "capacityfactor",
            indexed by gridcell x and y coordinates (e.g. "rlat" and "rlon").

    Returns:
        xarray.Dataset: The dataset contains two arrays:
            1. capacityfactors by gridcell
            2. eligible areas
        The dataset is trimmed to those gridcells with some eligble area.
    """
    area_df = pd.read_csv(path_to_eligible_areas)

    model_area_da = area_df.set_index(["rlat", "rlon", "GID_1"])["area"].to_xarray()



    ds = xr.open_dataset(path_to_turbine_cf)

    # The rlat/rlon cooridnates are always slightly different due to floating point errors,
    # so we have to re-interpolate the data to re-align the data
    ds_interp_to_area_coords = ds.interp(rlat=model_area_da.rlat, rlon=model_area_da.rlon, method="nearest")

    ds_interp_to_area_coords["area"] = model_area_da

    # Keep only those sites with eligible area > 0
    return (
        ds_interp_to_area_coords
        .stack(site=["rlat", "rlon"])
        .dropna(
            "site", subset=["area"], how="all"
        )
    )


def cluster_cf_data(timeseries_data, clustering_method, n_clusters):
    """
    Cluster timeseries data into n_clusters according to given clustering method

    Args:
        timeseries_data (numpy.ndarray or pandas.DataFrame):
            If clustering using SKLearn clustering algorithms: (n_samples, n_features).
            If clustering with tslearn TimeSeriesKMeans: (n_features, n_samples).
        clustering_method (str):
            One of "kmeans" (scikit-learn) , "spectral" (scikit-learn),
            or "timeseries_kmeans" (tslearn)
        n_clusters (int): number of clusters to create

    Returns:
        numpy.array: labels corresponding to clusters for each sample in n_samples.
    """
    if clustering_method == "kmeans":
        model = KMeans(
            n_clusters=n_clusters, **CLUSTERING_KWARGS["kmeans"]
        )
    elif clustering_method == "spectral":
        model = SpectralClustering(
            n_clusters=n_clusters, **CLUSTERING_KWARGS["spectral"]
        )
    elif clustering_method == "timeseries_kmeans":
        model = TimeSeriesKMeans(
            n_clusters=n_clusters, **CLUSTERING_KWARGS["timeseries_kmeans"]
        )

    return model.fit_predict(timeseries_data)


def get_region_cf_df(cf_ds, region_name=None):
    if region_name is not None:
        cf_ds = cf_ds.sel(GID_1=region_name).dropna("site", subset=["area"], how="all")
    cf_df = (
        cf_ds
        .capacityfactor
        .to_series()
        .droplevel("turbine")
        .unstack(["rlat", "rlon"])
    )

    # clean nan and tiny negative capacityfactors ready for clustering
    return cf_df.dropna(how="all").dropna(how="all", axis=1).clip(lower=0)


def get_region_area(cf_ds, region_name=None):
    if region_name is not None:
        cf_ds = cf_ds.sel(GID_1=region_name).dropna("site", subset=["area"], how="all")
    return (
        cf_ds
        .area
        .to_series()
    )


def get_weighted_average(timeseries_df, weights):
    return pd.Series(
        index=timeseries_df.index,
        data=np.average(timeseries_df, weights=weights.loc[timeseries_df.columns], axis=1)
    )


if __name__ == "__main__":
    create_clusters(
        path_to_turbine_cf=snakemake.input.turbine_cf,
        path_to_eligible_areas=snakemake.input.eligible_areas,
        clustering_method=snakemake.wildcards.clustering_method,
        model_region_config=snakemake.params.model_region_config,
        threads=snakemake.threads,
        path_to_timeseries_output=snakemake.output.timeseries_output,
        path_to_label_output=snakemake.output.label_output,
    )
