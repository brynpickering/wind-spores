def get_turbine_height_combos(heights):
    turbine_height_combos = []
    for height in heights:
        for turbine_name, turbine_config in config["study-turbines"].items():
            if (height >= turbine_config["height"]["min"]) & (height <= turbine_config["height"]["max"]):
                turbine_height_combos.append(f"{turbine_name}_{height}")
    return turbine_height_combos


def get_capacity_mw(wildcards):
    if wildcards.turbine_name in config["current-turbines"].keys():
        return config["current-turbines"][wildcards.turbine_name]["capacity_mw"]
    elif wildcards.turbine_name.rsplit("_", 1)[0] in config["study-turbines"].keys():
        return config["study-turbines"][wildcards.turbine_name.rsplit("_", 1)[0]]["capacity_mw"]


ruleorder: turbine_metrics_per_ch_unit > best_turbines_per_gridcell > area_weighted_aggregate_metrics

rule turbine_output:
    message: "Get {wildcards.turbine_name} turbine output for {wildcards.dataset_name} wind speeds"
    input:
        script = "scripts/vwf.py",
        wind_speed_coeffs = "build/{dataset_name}/A-Z-coeffs.nc",
        power_curves = config["data-sources"]["power-curves"]
    params:
        turbine_config = lambda wildcards: config["current-turbines"][wildcards.turbine_name],
        height = lambda wildcards: config["current-turbines"][wildcards.turbine_name]["height"]
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
        turbine_capacity = lambda wildcards: get_capacity_mw(wildcards),
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
        turbine_capacity = lambda wildcards: get_capacity_mw(wildcards),
        turbine_density = config["turbine-params"]["density"],
        dataset_config = lambda wildcards: config[wildcards.dataset_name],
        comparison_quantiles = [0.1, 0.25, 0.5, 0.75, 0.9, 1],
        comparison_periods = {"all": range(1, 13), "summer": (6, 7, 8), "winter": (1, 2, 12)},
        level = lambda wildcards: wildcards.level
    conda: "../envs/geo.yaml"
    wildcard_constraints:
        turbine_name = "|".join(set(config["current-turbines"].keys()) & set(get_turbine_height_combos(config["study-heights"]))),
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
    output: "build/{dataset_name}/top-current-turbines/{cf_or_mwh}-{ismasked}.nc"
    script: "../scripts/compare_turbines.py"


rule turbine_output_single:
    message: "Get {wildcards.turbine_name} turbine output for {wildcards.dataset_name} wind speeds"
    input:
        script = "scripts/vwf.py",
        wind_speed_coeffs = "build/{dataset_name}/A-Z-coeffs.nc",
        power_curves = config["data-sources"]["power-curves-single-turbine"]
    params:
        turbine_config = lambda wildcards: config["study-turbines"][wildcards.turbine_name.rsplit("_", 1)[0]],
        height = lambda wildcards: float(wildcards.turbine_name.rsplit("_", 1)[-1])
    conda: "../envs/default.yaml"
    output: "build/{dataset_name}/{turbine_name}.nc"
    wildcard_constraints:
        turbine_name = "|".join(get_turbine_height_combos(config["study-heights"])),
        dataset_name = "newa|cosmo-rea2"
    script: "../scripts/vwf.py"


rule best_study_turbines_per_gridcell:
    message:
        """
        Compare per-gridcell {wildcards.cf_or_mwh} performance of study turbines based
        on {wildcards.dataset_name} wind speeds over all dataset years,
        {wildcards.ismasked} to land eligibility
        """
    input:
        script = "scripts/compare_turbines.py",
        turbine_metrics = lambda wildcards: expand(
            "build/{{dataset_name}}/{turbine_name}-metrics-{{cf_or_mwh}}-{{ismasked}}.nc",
            turbine_name=get_turbine_height_combos(config["study-heights"]),
        )
    conda: "../envs/default.yaml"
    output: "build/{dataset_name}/top-study-turbines/{cf_or_mwh}-{ismasked}.nc"
    script: "../scripts/compare_turbines.py"
