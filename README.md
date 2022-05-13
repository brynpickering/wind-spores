# wind-spores

SFOE-funded project WindSPORES: Policy-relevant wind power deployment scenarios for Switzerland.

This repository contains the workflow to process data for the first half of the project: the analysis of wind data in Switzerland. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

The repository relies on the [`solar-and-wind-potentials`](https://github.com/calliope-project/solar-and-wind-potentials) submodule. Use this command to clone including the submodule:

    git clone --recurse-submodules git@github.com:brynpickering/wind-spores.git

Or, if you have already cloned the repository, first run:

    git submodule update --init --recursive

You need [conda](https://conda.io/docs/index.html) to run the analysis. Using conda, you can create a conda environment from within you can run it:

    conda env create -f environment.yaml

### Data to be retrieved manually
You will need some intermediate data to complete parts of the workflow, where access to this data via an online API is difficult or impossible. You will need to prepare this data yourself before initiating the analysis.

1. [COSMO-REA2 wind speed log law coefficients](https://reanalysis.meteo.uni-bonn.de/?Download_Data___COSMO-REA2). The coefficients "A" and "z" which describe the relationship between vertical height above ground and wind speed at every COSMO-REA2 gridcell, based on the equation ``speed = A log(height) - A log(z)``. For an example of the code to derive these coefficients, see the script `log_law_coefficients.py` where we undertake the proess for the NEWA dataset. The data should be split into one file per year and placed in `data/cosmo_rea2_A_Z/{year}.nc`.

2. Measured hourly turbine electricity generation data, for comparison with simulated data. This data is provided to us under a restricted license, so we cannot include it directly in this repository, nor can we host it elsewhere online. If you provide your own data, it should be a netcdf with turbine/farm site name as the data variables, and `datetime` (e.g. `YYYY-mm-DD HH:MM`) and `height` (hub height) as the data dimensions. To be placed in `data/turbine_power_output_update_for_windspores.nc`.

3. [EU-DEM slope data](https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1-0-and-derived-products/). "Full European Coverage" raster to be placed in `.data/eudem_slop_3035_europe.tif`.

Futhermore, there are two files needed by the `solar-and-wind-potentials` submodule (`European Settlement Map 2012, Release 2017, 100m` & `World Exclusive Economic Zones v10`) which should be placed in the data directory of the this workflow (i.e., not in the submodule data directory). Details on the data sources can be found in the submodule README. The capacity factor data referenced in that README is not necessary to run this workflow.

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
