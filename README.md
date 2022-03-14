# CALDERA Paper

## Organisation

This repository is organized as follow:

- [Dockerfile](./Dockerfile) is a docker instruction file that will install all requirements, both for CALDERA and for pre-processing and post-processing the data
- [Makefile](./Makefile) contains a list of make commands. `make build` builds to docker image using the `Dockerfile`. Further make commands are listed in the Analysis section
- [Raw](./Raw/): A folder where the raw data is
- [Data](./Data/): A folder where the pre-processed data is stored, before running CALDERA
- [Output](./Output/)`: The final output of each analysis
- [Analysis](./Analysis/): A folder which contains the scripts for the three analyses. More details in the Analysis section
- [Figures](./Figures/): A folder containing the Caldera algorithm, as well as our implementation of  COIN,  FACS and All Unitigs.
- [src](./src/): Scripts to generate the figures from the main paper and supplementary

## Running Caldera

Assuming we are on a docker container created via the [Dockerfile](./Dockerfile), we first create a test dataset.

```r
out=/Example/Output
raw=/Example/Raw
data=/Example/Data/
mkdir -p $out $raw/Genomes
R CMD BATCH Analysis/Example/generate.R ${out}/generate.out
```

This creates fasta files, one per genome, as well as a strains file. For more information, see `DBGWAS -h` to understand the proper format.

We then run the actual scripts:

- Run DBGWAS step1 (via the `-only1` parameter) to build the DBG.
- Run the `toMajor.py` script to prepare the format.
- Run the Caldera script with two threads (`-t 2`), on one community (`-C 1`).
- Run DBGWAS step 3 (`-skip1 -skip2`) to visualize the results. This creates an html file in Outpout/caldera/step3 where we can visualize the significant CCS.

```r
DBGWAS -strains ${raw}/strains -only1 -nb-cores 2 -output ${data} 
python3 src/CALDERA/bin/Pre-Process/toMajor.py -l ${data}
python3 src/CALDERA/bin/caldera-script -l ${data} -o ${out}/caldera/ -v -t 2 -C 1
DBGWAS -strains ${raw}/strains -skip1 -skip2  -nb-cores 2 \
    -output ${out}/caldera/ -nh 200 -caldera 
```

All those files can also be found [here](./Output/Example/).

## Using the docker container

Clone the repo with:

```sh
git clone https://github.com/HectorRDB/Caldera_ISMB.git
```

Then, build the dockerfile (you can use `make build`). To run the example, you can start an intereactive session via `make link`.
