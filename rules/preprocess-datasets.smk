import pandas as pd

localrules: download_newa_one_day  # this can't be parallelised due to throttling on the server-side


rule download_newa_one_day:
    message: "Download all NEWA meso-scale wind speeds for Switzerland on the date {wildcards.date}"
    output: temp("data/automatic/newa/{date}.nc")
    wildcard_constraints:
        date = "\d+-\d+-\d+"
    params:
        year = lambda wildcards: pd.to_datetime(wildcards.date).year
    conda: "../envs/shell.yaml"
    shell:  # wget snippet based on script from https://map.neweuropeanwindatlas.eu/
        """
        wget --content-on-error -O {output} http://opendap.neweuropeanwindatlas.eu:80/opendap/newa/NEWA_MESOSCALE_ATLAS/{params.year}/NEWA-{wildcards.date}.nc.nc?WS[0:1:47][0:1:4][396:1:481][524:1:647],WS10[0:1:47][396:1:481][524:1:647],time[0:1:47],height[0:1:4],west_east[524:1:647],south_north[396:1:481],XLON[396:1:481][524:1:647],XLAT[396:1:481][524:1:647],Times[0:1:47],crs
        ncks -O --mk_rec_dmn time {output} {output}
        """


rule newa_one_year:
    message: "Combine daily NEWA meso-scale wind speeds for Switzerland in the year {wildcards.year}"
    input:
        lambda wildcards: expand(
            "data/automatic/newa/{date}.nc",
            date=pd.date_range("{}-01-01".format(wildcards.year), "{}-12-31".format(wildcards.year), freq='D').astype(str)
        )
    output: protected("data/automatic/newa/{year}.nc")
    wildcard_constraints:
        year = "|".join(str(i) for i in range(2009, 2019))
    shadow: "minimal"
    conda: "../envs/shell.yaml"
    shell: "ncra --mro -O -d time,,,2,2 {input} {output}"


rule newa_annual_a_Z_interpolation:
    message: "Infer log-law coefficients from NEWA {wildcards.year} fixed height dataset, for input to the Renewables.ninja Virtual Wind Farm (VWF)"
    input:
        script = "scripts/log_law_coefficients.py",
        wind_data = rules.newa_one_year.output[0],
    params:
        height_weights = config["newa"]["height-weights"],
        x_name = config["newa"]["x-name"],
        y_name = config["newa"]["y-name"],
        lon_name = config["newa"]["lon-name"],
        lat_name = config["newa"]["lat-name"],
        wind_speed_var_name = "WS"
    output: temp("build/newa/{year}-A-Z-coeffs.nc")
    conda: "../envs/default.yaml"
    script: "../scripts/log_law_coefficients.py"


rule cosmo_switzerland_bbox:
    message: "Slice COSMO-REA2 dataset to extent of Switzerland (+ a buffer) for the year {wildcards.year}"
    input:
        script = "scripts/slice_cosmo.py",
        wind_speed = lambda wildcards: config["data-sources"]["cosmo-rea2-wind-speed-coeffs"].format(year=wildcards.year),
        cosmo_coords = config["data-sources"]["cosmo-rea2-coords"]
    params:
        bounding_box = config["ch-bounding-box"],
        cosmo_config = config["cosmo-rea2"]
    conda: "../envs/default.yaml"
    output: temp("build/cosmo-rea2/{year}-A-Z-coeffs.nc")
    script: "../scripts/slice_cosmo.py"


rule all_years_a_Z_coeffs:
    message: "Merge all years of {wildcards.dataset} log-law A & Z coefficients"
    input:
        lambda wildcards: expand(
            "build/{{dataset}}/{year}-A-Z-coeffs.nc",
            year=range(
                config[wildcards.dataset]["first-year"],
                config[wildcards.dataset]["final-year"] + 1
            )
        )
    output: "build/{dataset}/A-Z-coeffs.nc"
    wildcard_constraints:
        dataset = "((newa)|(cosmo-rea2))"
    run:
        import xarray as xr
        xr.concat(
            [xr.open_dataset(file) for file in input], dim="time"
        ).sortby("time").to_netcdf(output[0])


rule wind_speed_at_height:
    message: "Get wind speed at {wildcards.height}m above ground, using {wildcards.dataset} log-law coefficients"
    input:
        script = "scripts/coeffs_to_wind_speed.py",
        coeffs = "build/{dataset}/A-Z-coeffs.nc"
    conda: "../envs/default.yaml"
    output: "build/{dataset}/wind-speed-{height}m.nc"
    params:
        lon_name = lambda wildcards: config[wildcards.dataset]["lon-name"],
        lat_name = lambda wildcards: config[wildcards.dataset]["lat-name"]
    wildcard_constraints:
        dataset = "((newa)|(cosmo-rea2))"
    script: "../scripts/coeffs_to_wind_speed.py"


rule measured_cf:
    message: "Translate measured {wildcards.turbine_id} output at {wildcards.turbine_site} to capacity factor"
    input:
        script = "scripts/measured_turbine_output_to_cf.py",
        turbine_output = config["data-sources"]["turbine-output"],
        turbines_per_site = config["data-sources"]["turbines-per-site"],
        current_turbine_sites = config["data-sources"]["turbine-sites"]
    conda: "../envs/geo.yaml"
    params:
        turbine_config = lambda wildcards: config["current-turbines"][wildcards.turbine_id]
    output: "build/measured/{turbine_site}/{turbine_id}.nc"
    script: "../scripts/measured_turbine_output_to_cf.py"
