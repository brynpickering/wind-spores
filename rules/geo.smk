import xarray as xr

ruleorder: site_specific_data > regrid

rule cosmo_crs:
    message: "Get the geographic projection of COSMO-REA2 data"
    input:
        script = "scripts/crs.py",
    params:
        pole_latitude = config["cosmo-rea2"]["pole-latitude"],
        pole_longitude = config["cosmo-rea2"]["pole-longitude"],
        data_source = "cosmo-rea2"
    conda: "../envs/geo.yaml"
    output: "build/cosmo-rea2/crs.txt"
    script: "../scripts/crs.py"


rule newa_crs:
    message: "Get the geographic projection of NEWA data"
    input:
        script = "scripts/crs.py",
        data = "data/automatic/newa/2010.nc"
    params:
        data_source = "newa"
    conda: "../envs/geo.yaml"
    output: "build/newa/crs.txt"
    script: "../scripts/crs.py"


rule dataset_polygons:
    message: "Get polygons representing the points of {wildcards.dataset}, in the EPSG:4326 projection"
    input:
        script = "scripts/points_to_polys.py",
        data = "build/{dataset}/A-Z-coeffs.nc",
        crs = "build/{dataset}/crs.txt"
    params:
        config = lambda wildcards: config[wildcards.dataset],
        output_driver = "GeoJSON"
    conda: "../envs/geo.yaml"
    output: "build/{dataset}/polys.geojson"
    script: "../scripts/points_to_polys.py"


rule site_specific_data:
    message: "Interpolate gridded {wildcards.dataset} {wildcards.timeseries} data to coordinates of {wildcards.turbine_site} turbine site"
    input:
        script = "scripts/data_at_site.py",
        crs = lambda wildcards: "build/{}/crs.txt".format(wildcards.dataset.split("-gridded-to-")[-1]),
        modelled_data = "build/{dataset}/{timeseries}.nc",
        current_turbine_sites = config["data-sources"]["turbine-sites"]
    params:
        x_name = lambda wildcards: config[wildcards.dataset]["x-name"],
        y_name = lambda wildcards: config[wildcards.dataset]["y-name"],
        dataset_name = lambda wildcards: config[wildcards.dataset]["long-name"]
    conda: "../envs/geo.yaml"
    output: "build/{dataset}/{turbine_site}/{timeseries}.nc"
    wildcard_constraints:
        turbine_site = "((martigny)|(st_brais))",
        dataset = "((newa)|(cosmo-rea2)|(cosmo-rea2-gridded-to-newa))"
    script: "../scripts/data_at_site.py"


rule all_site_specific_data:
    message: "Combine all {wildcards.timeseries} datasets interpolated to coordinates of {wildcards.turbine_site} turbine site"
    input:
        expand(
            "build/{dataset}/{{turbine_site}}/{{timeseries}}.nc",
            dataset=["newa", "cosmo-rea2", "cosmo-rea2-gridded-to-newa", "measured"]
        ),
    output: "build/{turbine_site}/{timeseries}.nc"
    wildcard_constraints:
        turbine_site = "((martigny)|(st_brais))"
    run:
        import xarray as xr
        xr.concat(
            [xr.open_dataset(file) for file in input], dim="dataset"
        ).to_netcdf(output[0])


rule regrid:
    message: "regrid {wildcards.timeseries} {wildcards.in_dataset_name} to that of {wildcards.out_dataset_name}"
    input:
        script = "scripts/regrid.py",
        in_data = "build/{in_dataset_name}/{timeseries}.nc",
        out_data = "build/{out_dataset_name}/{timeseries}.nc"
    params:
        in_config = lambda wildcards: config[wildcards.in_dataset_name],
        out_config = lambda wildcards: config[wildcards.out_dataset_name]
    wildcard_constraints:
        in_dataset_name = "((newa)|(cosmo-rea2))",
        out_dataset_name = "((newa)|(cosmo-rea2))",
    conda: "../envs/geo.yaml"
    output: "build/{in_dataset_name}-gridded-to-{out_dataset_name}/{timeseries}.nc"
    script: "../scripts/regrid.py"


rule ch_shape_zip:
    message: "Download Swiss border shape as ZIP file"
    params:
        url = config["data-sources"]["ch-shape"]
    conda: "../envs/shell.yaml"
    output: "data/automatic/che.gpkg.zip"
    shell: "curl -sLo {output} '{params.url}'"


rule eligible_area_per_gridcell_and_ch_level:
    message:
        """
        Assign eligible wind areas to each {wildcards.dataset_name} grid and
        Swiss level {wildcards.level} administrative unit
        """
    input:
        script = "scripts/eligible_area_per_gridcell_and_region.py",
        polygons = "build/{dataset_name}/polys.geojson",
        eligible_land = "build/technically-eligible-land-wind.tif",
        ch_shape = rules.ch_shape_zip.output[0]
    params:
        dataset_config = lambda wildcards: config[wildcards.dataset_name]
    wildcard_constraints:
        level = "0|1|2|3",
        dataset_name = "((newa)|(cosmo-rea2))"
    conda: "../envs/geo.yaml"
    output: "build/{dataset_name}-{cf_or_mwh}-gridded-to-ch-level-{level}/eligible-areas.csv"
    script: "../scripts/eligible_area_per_gridcell_and_region.py"


rule eligible_area_per_gridcell_and_model_region:
    message:
        """
        Assign eligible wind areas to each {wildcards.dataset_name} grid and
        Euro-Calliope swiss region
        """
    input:
        script = "scripts/eligible_area_per_gridcell_and_region.py",
        eligible_areas = "build/{dataset_name}-{cf_or_mwh}-gridded-to-ch-level-1/eligible-areas.csv",
        model_region_mapping = config["data-sources"]["model-region-mapping"]
    wildcard_constraints:
        dataset_name = "((newa)|(cosmo-rea2))"
    conda: "../envs/geo.yaml"
    output: "build/{dataset_name}-{cf_or_mwh}-gridded-to-ch-model-regions/eligible-areas.csv"
    script: "../scripts/eligible_area_per_gridcell_and_region.py"


rule area_weighted_aggregate_metrics:
    message:
        """
        Aggregate {wildcards.dataset_name} gridded {wildcards.timeseries} {wildcards.cf_or_mwh} data
        to Swiss level {wildcards.level} administrative units,
        using a technically available area weighted average
        """
    input:
        script = "scripts/aggregate_to_ch_units.py",
        timeseries_data = "build/{dataset_name}/{timeseries}.nc",
        eligible_areas = rules.eligible_area_per_gridcell_and_ch_level.output[0],
    params:
        dataset_config = lambda wildcards: config[wildcards.dataset_name]
    wildcard_constraints:
        level = "0|1|2|3",
        cf_or_mwh = "CF|MWh",
        dataset_name = "((newa)|(cosmo-rea2))",
        timeseries = "[^/]*"
    conda: "../envs/default.yaml"
    output: "build/{dataset_name}-{cf_or_mwh}-gridded-to-ch-level-{level}/{timeseries}.nc"
    script: "../scripts/aggregate_to_ch_units.py"


rule cluster_model_regions:
    message: "Cluster {wildcards.dataset_name} gridded {wildcards.timeseries} CF within model regions using {wildcards.clustering_method}"
    input:
        script = "scripts/cluster_cf_to_model_regions.py",
        turbine_cf = "build/{dataset_name}/{timeseries}.nc",
        eligible_areas = "build/{dataset_name}-CF-gridded-to-ch-model-regions/eligible-areas.csv",
    params:
        model_region_config = config["model-regions"]
    wildcard_constraints:
        dataset_name = "((newa)|(cosmo-rea2))",
    conda: "../envs/default.yaml"
    threads: 5
    output:
        timeseries_output = "build/{dataset_name}-CF-gridded-to-ch-model-regions/{timeseries}-clustered-CF-{clustering_method}.csv",
        label_output = "build/{dataset_name}-CF-gridded-to-ch-model-regions/{timeseries}-cluster-labels-{clustering_method}.csv"
    script: "../scripts/cluster_cf_to_model_regions.py"


rule test_clusters:
    message: """
        Test different number of {wildcards.clustering_method} clustered regions across CH,
        based on {wildcards.turbine_name} capacity factors derived from {wildcards.dataset_name} weather data
        """
    input:
        script = "scripts/test_clustering_options.py",
        turbine_cf = "build/{dataset_name}/{turbine_name}.nc",
        eligible_areas = "build/{dataset_name}-CF-gridded-to-ch-model-regions/eligible-areas.csv",
    wildcard_constraints:
        dataset_name = "((newa)|(cosmo-rea2))",
    conda: "../envs/default.yaml"
    threads: 5
    output: "build/{dataset_name}/{turbine_name}-cluster-labels-{clustering_method}.csv"
    script: "../scripts/test_clustering_options.py"
