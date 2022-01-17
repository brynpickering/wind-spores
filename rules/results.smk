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
        top_turbines = "build/{dataset}/top-turbines-CF-unmasked.nc",
        ch_shape = rules.ch_shape_zip.output[0],
        polys = "build/{dataset}/polys.geojson",
    params:
        plot_crs = "EPSG:3035",
        dataset_config = lambda wildcards: config[wildcards.dataset],
        turbine_config = config["turbines"],
        level = None
    conda: "../envs/vis.yaml"
    wildcard_constraints:
        suffix = "((png)|(pdf))"
    output: "build/results/compare_top_turbines_per_{dataset}_gridcell.{suffix}"
    script: "../scripts/compare_best_turbines_per_gridcell.py"


rule plot_turbine_comparison_per_region:
    message: "Plot top turbines per Swiss level {wildcards.level} administrative units, based on {wildcards.cf_or_mwh} derived from {wildcards.ismasked} {wildcards.dataset_name} weather data"
    input:
        script = "scripts/compare_best_turbines_per_gridcell.py",
        top_turbines = "build/{dataset_name}-CF-gridded-to-ch-level-{level}/top-turbines-{cf_or_mwh}-{ismasked}.nc",
        ch_shape = rules.ch_shape_zip.output[0],
    params:
        plot_crs = "EPSG:3035",
        dataset_config = lambda wildcards: config[wildcards.dataset_name],
        turbine_config = config["turbines"],
        level = lambda wildcards: wildcards.level
    conda: "../envs/vis.yaml"
    wildcard_constraints:
        suffix = "((png)|(pdf))"
    output: "build/results/compare_top_turbines_per_CH_level{level}_from_{dataset_name}-{cf_or_mwh}-{ismasked}.{suffix}"
    script: "../scripts/compare_best_turbines_per_gridcell.py"
