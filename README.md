# wind-spores

SFOE-funded project WindSPORES: Policy-relevant wind power deployment scenarios for Switzerland.

This repository contains the workflow to process data for the first half of the project: the analysis of wind data in Switzerland. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

You need [conda](https://conda.io/docs/index.html) to run the analysis. Using conda, you can create a conda environment from within you can run it:

    conda env create -f environment.yaml

You also will need some intermediate data to complete parts of the workflow, where access to this data via an online API is difficult or impossible. You will need to prepare this data yourself before initiating the analysis, then point to the data filepaths in `config/full.yaml`:

1. COSMO-REA2 wind speed log law coefficients (`cosmo-rea2-wind-speed-coeffs`). The coefficients "A" and "z" which describe the relationship between vertical height above ground and wind speed at every COSMO-REA2 gridcell, based on the equation ``speed = A log(height) - A log(z)``. For an example of the code to derive these coefficients, see the script `log_law_coefficients.py` where we undertake the proess for the NEWA dataset. The data should be split into one file per year.

2. Measured hourly turbine electricity generation data, for comparison with simulated data (`turbine-output`). This data is provided to us under a restricted license, so we cannot include it directly in this repository, nor can we host it elsewhere online. If you provide your own data, it should be a netcdf with turbine/farm site name as the data variables, and `datetime` (e.g. `YYYY-mm-DD HH:MM`) and `height` (hub height) as the data dimensions.

## Run the analysis

    snakemake --configfile config/full.yaml --use-conda

This will run all analysis steps to reproduce results and figures used in reporting.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps, and if you have `dot` installed, run:

    snakemake --rulegraph | dot -Tpdf > dag.pdf



## Be notified of build successes or fails

  As the execution of this workflow may take a while, you can be notified whenever the execution terminates either successfully or unsuccessfully. Notifications are sent by email. To activate notifications, add the email address of the recipient to the configuration key `email`. You can add the key to your configuration file, or you can run the workflow the following way to receive notifications:

      snakemake --use-conda --config email=<your-email>

## Repo structure

* `scripts`: contains the Python source code as scripts
* `rules`: contains Snakemake rule definitions
* `envs`: contains execution environments
* `config`: configurations used in the study
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)

## License

The code in this repo is MIT licensed, see `./LICENSE.md`.
