rule cosmo_dataset_polygons:
    message: "Get polygons representing the points of cosmo-rea2, in the projection of the dataset"
    input:
        script = "scripts/points_to_polys.py",
        data = "build/cosmo-rea2/2010-A-Z-coeffs.nc"
    params:
        data_source = "cosmo-rea2",
        config = config["cosmo-rea2"],
        output_driver = "GeoJSON"
    conda: "../envs/geo.yaml"
    output: "build/cosmo-rea2/polys.geojson"
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
    message: "regrid {wildcards.dataset} {wildcards.in_dataset_name} to that of {wildcards.out_dataset_name}"
    input:
        script = "scripts/regrid.py",
        in_data = "build/{in_dataset_name}/{dataset}",
        out_data = "build/{out_dataset_name}/{dataset}"
    params:
        in_config = lambda wildcards: config[wildcards.in_dataset_name],
        out_config = lambda wildcards: config[wildcards.out_dataset_name]
    conda: "../envs/geo.yaml"
    output: "build/regridded/{in_dataset_name}-gridded-to-{out_dataset_name}/{dataset}"
    script: "../scripts/regrid.py"
