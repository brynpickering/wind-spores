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


rule newa_dataset_polygons:
    message: "Get polygons representing the points of NEWA, in the projection of the dataset"
    input:
        script = "scripts/points_to_polys.py",
        data = "data/automatic/newa/2010.nc"
    params:
        data_source = "newa",
        config = config["newa"],
        output_driver = "GeoJSON"
    conda: "../envs/geo.yaml"
    output: "build/newa/polys.geojson"
    script: "../scripts/points_to_polys.py"


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
