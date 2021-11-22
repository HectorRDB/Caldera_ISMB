# CALDERA Paper

This repository is organized as follow:

- `Dockerfile` is a docker instruction file that will install all requirements, both for CALDERA and for pre-processing and post-processing the data
- `Makefile` contains a list of make commands. `make build` builds to docker image using the `Dockerfile`. Further make commands are listed in the Analysis section
- `Raw`: A folder where the raw data is
- `Data`: A folder where the pre-processed data is stored, before running CALDERA
- `Output`: The final output of each analysis
- `Analysis`: A folder which contains the scripts for the three analyses. More details in the Analysis section
- `scr`: A folder containing the Caldera algorithm, as well as our implementation of  COIN,  FACS and All Unitigs.
- `Figures`: Scripts to generate the figures from the main paper and supplementary

You can use the Analysis/Example folder to run a simple example that should run under 5 minutes. The output from that run is already present in Output/Example.
