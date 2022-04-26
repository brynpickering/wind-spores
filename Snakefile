configfile: "config/default.yaml"

include: "./rules/preprocess-datasets.smk"
include: "./rules/sync.smk"
include: "./rules/geo.smk"
include: "./rules/turbines.smk"
include: "./rules/results.smk"
include: "./rules/eligible-land.smk"

onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'wind-spores succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'wind-spores failed' {config[email]}")


rule all:
    message: "Run entire analysis and produce report figures."
    input:
        "build/results/compare_cosmo_newa_annual_average_100m.png",
        "build/results/compare_specific_sites_cosmo_newa.png",
        "build/results/compare_top_study_turbines_per_cosmo-rea2_gridcell.png",
        "build/results/compare_study_turbines_swiss_average_from_cosmo-rea2-CF.png",
        "build/results/compare-kmeans-clustering-levels-cosmo-rea2-vestas_v110_2000_120-CF.png",
        "build/results/compare-kmeans-clustering-cosmo-rea2-CF.png",
        "build/results/compare-spectral-clustering-levels-cosmo-rea2-vestas_v110_2000_120-CF.png",


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """


rule test:
    conda: "envs/test.yaml"
    output: "build/test-report.html"
    shell:
        "py.test --html={output} --self-contained-html"
