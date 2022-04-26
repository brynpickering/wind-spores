rule plot_compare_cosmo_newa_annual_average:
    message: "Create {wildcards.suffix} figure to compare wind speeds at {wildcards.height}m above ground, according to COSMO-REA2 and NEWA"
    input:
        script = "scripts/compare_cosmo_newa_annual_average.py",
        cosmo_wind_speed = "build/cosmo-rea2/wind-speed-{height}m.nc",
        newa_wind_speed = "build/newa/wind-speed-{height}m.nc",
        cosmo_regridded_wind_speed = "build/cosmo-rea2-gridded-to-newa/wind-speed-{height}m.nc",
        cosmo_polys = "build/cosmo-rea2/polys.geojson",
        newa_polys = "build/newa/polys.geojson",
        ch_shape = rules.ch_shape_zip.output[0]
    params:
        plot_crs = "EPSG:3035",
        config = config,
    conda: "../envs/vis.yaml"
    output: "build/results/compare_cosmo_newa_annual_average_{height}m.{suffix}"
    wildcard_constraints:
        suffix = "((png)|(pdf))"
    script: "../scripts/compare_cosmo_newa_annual_average.py"


rule plot_compare_specific_sites_cosmo_newa:
    message: "Create {wildcards.suffix} figure to compare median capacity factors between measured and simulated datasets"
    input:
        script = "scripts/compare_specific_sites_cosmo_newa.py",
        martigny_data = "build/martigny/enercon_e82_2000_99.nc",
        st_brais_data = "build/st_brais/enercon_e82_2000_78.nc"
    conda: "../envs/vis.yaml"
    output: "build/results/compare_specific_sites_cosmo_newa.{suffix}"
    wildcard_constraints:
        suffix = "((png)|(pdf))"
    script: "../scripts/compare_specific_sites_cosmo_newa.py"


rule plot_turbine_comparison_per_gridcell:
    message: "Plot top turbines per {wildcards.dataset} gridcell"
    input:
        script = "scripts/compare_best_turbines_per_gridcell.py",
        top_turbines = "build/{dataset}/top-{turbine_group}-turbines/CF-unmasked.nc",
        ch_shape = rules.ch_shape_zip.output[0],
        polys = "build/{dataset}/polys.geojson",
    params:
        plot_crs = "EPSG:3035",
        dataset_config = lambda wildcards: config[wildcards.dataset],
        turbine_config = lambda wildcards: config[f"{wildcards.turbine_group}-turbines"],
        level = None
    conda: "../envs/vis.yaml"
    wildcard_constraints:
        suffix = "((png)|(pdf))"
    output: "build/results/compare_top_{turbine_group}_turbines_per_{dataset}_gridcell.{suffix}"
    script: "../scripts/compare_best_turbines_per_gridcell.py"


rule plot_turbine_comparison_per_region:
    message: "Plot top turbines per Swiss level {wildcards.level} administrative units, based on {wildcards.cf_or_mwh} derived from {wildcards.ismasked} {wildcards.dataset_name} weather data"
    input:
        script = "scripts/compare_best_turbines_per_gridcell.py",
        top_turbines = "build/{dataset_name}-CF-gridded-to-ch-level-{level}/top-{turbine_group}-turbines/{cf_or_mwh}-{ismasked}.nc",
        ch_shape = rules.ch_shape_zip.output[0],
    params:
        plot_crs = "EPSG:3035",
        dataset_config = lambda wildcards: config[wildcards.dataset_name],
        turbine_config = lambda wildcards: config[f"{wildcards.turbine_group}-turbines"],
        level = lambda wildcards: wildcards.level
    conda: "../envs/vis.yaml"
    wildcard_constraints:
        suffix = "((png)|(pdf))"
    output: "build/results/compare_top_{turbine_group}_turbines_per_CH_level{level}_from_{dataset_name}-{cf_or_mwh}-{ismasked}.{suffix}"
    script: "../scripts/compare_best_turbines_per_gridcell.py"


rule plot_annual_data_study_turbines:
    message: "Plot swiss averaged data for study turbines, based on Capacity factors derived from {wildcards.dataset_name} weather data"
    input:
        script = "scripts/compare_study_turbines.py",
        cfs = expand(
            "build/{{dataset_name}}-CF-gridded-to-ch-level-0/{turbine_name}.nc",
            turbine_name=get_turbine_height_combos(config["study-heights"])
        ),
        power_curves = config["data-sources"]["power-curves-single-turbine"]
    params:
        turbine_config = config["study-turbines"]
    conda: "../envs/vis.yaml"
    wildcard_constraints:
        suffix = "((png)|(pdf))"
    output: "build/results/compare_study_turbines_swiss_average_from_{dataset_name}-CF.{suffix}"
    script: "../scripts/compare_study_turbines.py"

