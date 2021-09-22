rule ch_shape_zip:
    message: "Download Swiss border shape as ZIP file"
    params:
        url = config["data-sources"]["ch-shape"]
    conda: "../envs/shell.yaml"
    output: "data/automatic/che.gpkg.zip"
    shell: "curl -sLo {output} '{params.url}'"


rule plot_compare_cosmo_newa_annual_average:
    message: "Create {wildcards.suffix} figure to compare wind speeds at {wildcards.height}m above ground, according to COSMO-REA2 and NEWA"
    input:
        script = "scripts/compare_cosmo_newa_annual_average.py",
        cosmo_wind_speed = "build/cosmo-rea2/wind-speed-{height}m.nc",
        newa_wind_speed = "build/newa/wind-speed-{height}m.nc",
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

