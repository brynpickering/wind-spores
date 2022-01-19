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


rule area_weighted_aggregate_metrics:
    message:
        """
        Aggreagate {wildcards.dataset_name} gridded {wildcards.timeseries} {wildcards.cf_or_mwh} data
        to Swiss level {wildcards.level} administrative units,
        using a technically available area weighted average
        """
    input:
        script = "scripts/aggregate_to_ch_units.py",
        timeseries_data = "build/{dataset_name}/{timeseries}.nc",
        polygons = "build/{dataset_name}/polys.geojson",
        eligible_land = "build/technically-eligible-land.tif",
        ch_shape = rules.ch_shape_zip.output[0]
    params:
        dataset_config = lambda wildcards: config[wildcards.dataset_name]
    wildcard_constraints:
        level = "0|1|2|3",
        cf_or_mwh = "CF|MWh",
        dataset_name = "((newa)|(cosmo-rea2))",
    conda: "../envs/geo.yaml"
    output: "build/{dataset_name}-{cf_or_mwh}-gridded-to-ch-level-{level}/{timeseries}.nc"
    script: "../scripts/aggregate_to_ch_units.py"
