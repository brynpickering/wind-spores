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
