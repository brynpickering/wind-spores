configfile: "config/renewables-potentials.yaml"

module renewablespotentials:
    snakefile: "../solar-and-wind-potentials/rules/data-preprocessing.smk"
    config: config["renewables-potentials"]

use rule * from renewablespotentials


rule category_of_technical_eligibility:
    message:
        "Determine upper bound surface eligibility for renewables based on land cover, slope, bathymetry, and settlements."
    input:
        src = "scripts/technical_eligibility.py",
        land_cover = rules.land_cover_in_europe.output[0],
        slope_pv = "build/slope-europe-pv.tif",
        slope_wind = "build/slope-europe-wind.tif",
        building_share = rules.settlements.output.buildings,
        urban_green_share = rules.settlements.output.urban_greens
    params:
        max_building_share = config["renewables-potentials"]["parameters"]["max-building-share"],
        max_urban_green_share = config["renewables-potentials"]["parameters"]["max-urban-green-share"],
        slope_threshold = config["renewables-potentials"]["parameters"]["max-slope-pixel-fraction-threshold"]
    output:
        "build/technically-eligible-land.tif"
    conda: "../envs/geo.yaml"
    script: "../scripts/technical_eligibility.py"


rule tech_slope_thresholds:
    message: "Create binary raster for {wildcards.tech}, whose land use is limited by slope using {threads} threads"
    input:
        scripts = "scripts/slopes.py",
        slopes_in_europe = config["renewables-potentials"]["data-sources"]["slope"]
    output: temp("build/data/eudem_slop_3035_europe_{tech}.tif")
    threads: config["renewables-potentials"]["snakemake"]["max-threads"]
    params:
        max_slope = lambda wildcards: config["renewables-potentials"]["parameters"]["max-slope"][wildcards.tech],
        max_threads = config["renewables-potentials"]["snakemake"]["max-threads"]
    conda: "../envs/geo.yaml"
    script: "../scripts/slopes.py"


rule slope_thresholds_warped_to_land_cover:
    message: "Warp {wildcards.tech} land availability according to slope to resolution of study using {threads} threads."
    input:
        land_cover = rules.land_cover_in_europe.output,
        slope_threshold = "build/data/eudem_slop_3035_europe_{tech}.tif"
    output: "build/slope-europe-{tech}.tif"
    threads: config["renewables-potentials"]["snakemake"]["max-threads"]
    conda: "../envs/geo.yaml"
    shell:
        """
        rio warp {input.slope_threshold} -o {output} --like {input.land_cover} \
        --resampling average --threads {threads}
        """
