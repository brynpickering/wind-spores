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


def get_config(dataset_name, config, config_item):
    if dataset_name in config.keys():
        return config[dataset_name][config_item]
    else:
        if config_item == "long-name":
            new_dataset_name = dataset_name.split("-gridded-to-")[0]
            return f"Regridded {config[new_dataset_name][config_item]}"
        else:
            new_dataset_name = dataset_name.split("-gridded-to-")[-1]
            return config[new_dataset_name][config_item]


rule site_specific_data:
    message: "Interpolate gridded {wildcards.dataset} {wildcards.timeseries} data to coordinates of {wildcards.turbine_site} turbine site"
    input:
        script = "scripts/data_at_site.py",
        crs = lambda wildcards: "build/{}/crs.txt".format(wildcards.dataset.split("-gridded-to-")[-1]),
        modelled_data = "build/{dataset}/{timeseries}.nc",
        current_turbine_sites = config["data-sources"]["turbine-sites"]
    params:
        x_name = lambda wildcards: get_config(wildcards.dataset, config, "x-name"),
        y_name = lambda wildcards: get_config(wildcards.dataset, config, "y-name"),
        dataset_name = lambda wildcards: get_config(wildcards.dataset, config, "long-name")
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
