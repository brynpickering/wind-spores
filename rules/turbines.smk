ruleorder: turbine_metrics_per_ch_unit > best_turbines_per_gridcell > area_weighted_aggregate_metrics

rule turbine_output:
    message: "Get {wildcards.turbine_name} turbine output for {wildcards.dataset_name} wind speeds"
    input:
        script = "scripts/vwf.py",
        wind_speed_coeffs = "build/{dataset_name}/A-Z-coeffs.nc",
        power_curves = config["data-sources"]["power-curves"]
    params:
        turbine_config = lambda wildcards: config["current-turbines"][wildcards.turbine_name]
    conda: "../envs/default.yaml"
    output: "build/{dataset_name}/{turbine_name}.nc"
    wildcard_constraints:
        turbine_name = "|".join(config["current-turbines"].keys()),
        dataset_name = "((newa)|(cosmo-rea2))"
    script: "../scripts/vwf.py"


rule turbine_metrics_per_gridcell:
    message:
        """
        Compare per-gridcell {wildcards.cf_or_mwh} performance
        of {wildcards.turbine_name} turbine based on {wildcards.dataset_name}
        wind speeds over all dataset years, {wildcards.ismasked} to land eligibility
        """
    input:
        script = "scripts/turbine_metrics.py",
        turbine_output = rules.turbine_output.output[0],
        polygons = "build/{dataset_name}/polys.geojson",
        eligible_land = lambda wildcards: "build/technically-eligible-land.tif" if wildcards.ismasked == "masked" else []
    params:
        turbine_capacity = lambda wildcards: config["current-turbines"][wildcards.turbine_name]["capacity_mw"],
        turbine_density = config["turbine-params"]["density"],
        dataset_config = lambda wildcards: config[wildcards.dataset_name],
        comparison_quantiles = [0.1, 0.25, 0.5, 0.75, 0.9, 1],
        comparison_periods = {"all": range(1, 13), "summer": (6, 7, 8), "winter": (1, 2, 12)},
        level = None
    conda: "../envs/geo.yaml"
    wildcard_constraints:
        ismasked = "masked|unmasked",
        cf_or_mwh = "CF|MWh"
    output: "build/{dataset_name}/{turbine_name}-metrics-{cf_or_mwh}-{ismasked}.nc"
    script: "../scripts/turbine_metrics.py"


rule turbine_metrics_per_ch_unit:
    message:
        """
        Compare level {wildcards.level} administrative unit {wildcards.cf_or_mwh} performance
        of {wildcards.turbine_name} turbine based on {wildcards.dataset_name}
        wind speeds over all dataset years, {wildcards.ismasked} to land eligibility
        """
    input:
        script = "scripts/turbine_metrics.py",
        turbine_output = "build/{dataset_name}-CF-gridded-to-ch-level-{level}/{turbine_name}.nc",
        polygons = rules.ch_shape_zip.output[0],
        eligible_land = lambda wildcards: "build/technically-eligible-land.tif" if wildcards.ismasked == "masked" else []
    params:
        turbine_capacity = lambda wildcards: config["current-turbines"][wildcards.turbine_name]["capacity_mw"],
        turbine_density = config["turbine-params"]["density"],
        dataset_config = lambda wildcards: config[wildcards.dataset_name],
        comparison_quantiles = [0.1, 0.25, 0.5, 0.75, 0.9, 1],
        comparison_periods = {"all": range(1, 13), "summer": (6, 7, 8), "winter": (1, 2, 12)},
        level = lambda wildcards: wildcards.level
    conda: "../envs/geo.yaml"
    wildcard_constraints:
        turbine_name = "|".join(config["current-turbines"].keys()),
        ismasked = "masked|unmasked",
        cf_or_mwh = "CF|MWh"
    output: "build/{dataset_name}-CF-gridded-to-ch-level-{level}/{turbine_name}-metrics-{cf_or_mwh}-{ismasked}.nc"
    script: "../scripts/turbine_metrics.py"


rule best_turbines_per_gridcell:
    message:
        """
        Compare per-gridcell {wildcards.cf_or_mwh} performance of turbines based
        on {wildcards.dataset_name} wind speeds over all dataset years,
        {wildcards.ismasked} to land eligibility
        """
    input:
        script = "scripts/compare_turbines.py",
        turbine_metrics = lambda wildcards: expand(
            "build/{{dataset_name}}/{turbine_name}-metrics-{{cf_or_mwh}}-{{ismasked}}.nc",
            turbine_name=config["current-turbines"].keys(),
        )
    conda: "../envs/default.yaml"
    output: "build/{dataset_name}/top-turbines-{cf_or_mwh}-{ismasked}.nc"
    script: "../scripts/compare_turbines.py"
