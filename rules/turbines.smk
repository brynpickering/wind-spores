rule turbine_output:
    message: "Get {wildcards.turbine_name} turbine output for {wildcards.dataset_name} wind speeds"
    input:
        script = "scripts/vwf.py",
        wind_speed_coeffs = "build/{dataset_name}/A-Z-coeffs.nc",
        power_curves = config["data-sources"]["power-curves"]
    params:
        turbine_config = lambda wildcards: config["turbines"][wildcards.turbine_name]
    conda: "../envs/default.yaml"
    output: "build/{dataset_name}/{turbine_name}.nc"
    wildcard_constraints:
        turbine_name = "|".join(config["turbines"].keys()),
        dataset_name = "((newa)|(cosmo-rea2))"
    script: "../scripts/vwf.py"


rule turbine_metrics_per_gridcell:
    message: "Compare per-gridcell performance of {wildcards.turbine_name} turbine based on {wildcards.dataset_name} wind speeds over all dataset years"
    input:
        script = "scripts/turbine_metrics.py",
        turbine_output = rules.turbine_output.output[0]
    params:
        comparison_quantiles = [0.1, 0.25, 0.5, 0.75, 0.9, 1],
        comparison_periods = {"all": range(1, 13), "summer": (6, 7, 8), "winter": (1, 2, 12)}
    conda: "../envs/default.yaml"
    output: "build/{dataset_name}/{turbine_name}-metrics.nc"
    script: "../scripts/turbine_metrics.py"


rule best_turbines_per_gridcell:
    message: "Compare per-gridcell performance of turbines based on {wildcards.dataset_name} wind speeds over all dataset years"
    input:
        script = "scripts/compare_turbines.py",
        turbine_metrics = lambda wildcards: expand(
            "build/{{dataset_name}}/{turbine_name}-metrics.nc",
            turbine_name=config["turbines"].keys(),
        )
    conda: "../envs/default.yaml"
    output: "build/{dataset_name}/top-turbines.nc"
    script: "../scripts/compare_turbines.py"
