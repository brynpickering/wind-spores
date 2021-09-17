PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --citeproc"

configfile: "config/default.yaml"

include: "./rules/preprocess-datasets.smk"
include: "./rules/turbines.smk"


onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'wind-spores succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'wind-spores failed' {config[email]}")

rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/report.html",
        "build/test-report.html"


rule run:
    message: "Runs the demo model."
    input: "scripts/model.py"
    params:
        slope = config["slope"],
        x0 = config["x0"]
    output: "build/results.pickle"
    conda: "envs/default.yaml"
    script: "scripts/model.py"


rule plot:
    message: "Visualises the demo results."
    input:
        scripts = "scripts/vis.py",
        results = rules.run.output
    output: "build/plot.png"
    conda: "envs/default.yaml"
    script: "scripts/vis.py"


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--self-contained --to html5"
    elif suffix == "pdf":
        return "--pdf-engine weasyprint"
    elif suffix == "docx":
        return []
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")


rule report:
    message: "Compile report.{wildcards.suffix}."
    input:
        "report/literature.yaml",
        "report/report.md",
        "report/pandoc-metadata.yaml",
        "report/apa.csl",
        "report/reset.css",
        "report/report.css",
        rules.plot.output
    params: options = pandoc_options
    output: "build/report.{suffix}"
    wildcard_constraints:
        suffix = "((html)|(pdf)|(docx))"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} report.md  --metadata-file=pandoc-metadata.yaml {params.options} \
        -o ../build/report.{wildcards.suffix}
        """


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
