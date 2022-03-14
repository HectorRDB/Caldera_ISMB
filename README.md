# CALDERA Paper

* [Organisation](#organisation)
* [Running Caldera](#running-caldera)
* [Using the docker container](#using-the-docker-container)
* [All options](#all-options)

## Organisation

This repository is organized as follow:

* [Dockerfile](./Dockerfile) is a docker instruction file that will install all requirements, both for CALDERA and for pre-processing and post-processing the data
* [Makefile](./Makefile) contains a list of make commands. `make build` builds to docker image using the `Dockerfile`. Further make commands are listed in the Analysis section
* [Raw](./Raw/): A folder where the raw data is
* [Data](./Data/): A folder where the pre-processed data is stored, before running CALDERA
* [Output](./Output/)`: The final output of each analysis
* [Analysis](./Analysis/): A folder which contains the scripts for the three analyses. More details in the Analysis section
* [Figures](./Figures/): A folder containing the Caldera algorithm, as well as our implementation of  COIN,  FACS and All Unitigs.
* [src](./src/): Scripts to generate the figures from the main paper and supplementary

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

* Run DBGWAS step1 (via the `-only1` parameter) to build the DBG.
* Run the `toMajor.py` script to prepare the format.
* Run the Caldera script with two threads (`-t 2`), on one community (`-C 1`).
* Run DBGWAS step 3 (`-skip1 -skip2`) to visualize the results. This creates an html file in Outpout/caldera/step3 where we can visualize the significant CCS.

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

## All options

```sh
usage: caldera-script [-h] -l LOC [-o OUTPUT] [-t THREADS] [-C COMMUNITIES]
                      [-P COMFILE] [--Lmax LMAX] [-s SMAX] [-b BATCH_SIZE]
                      [--alpha ALPHA] [-v] [--kmax KMAX] [--dfs] [--save_int]
                      [--restart]

Run the CALDERA script.

optional arguments:
  -h, --help            show this help message and exit
  -l LOC, --loc LOC     where to find the step 1 folder from DGWAS.
  -o OUTPUT, --output OUTPUT
                        Where to store the output. Default to loc/step2/
  -t THREADS, --threads THREADS
                        How many threads to use. Default to 1
  -C COMMUNITIES, --communities COMMUNITIES
                        Number of communities to find. Default to 3
  -P COMFILE, --comFile COMFILE
                        Location of the community assignement file. If none is
                        specified, default to k-means with k = communities
  --Lmax LMAX           Maximum size of each subgraph in bp. Default to 100
  -s SMAX, --Smax SMAX  Maximum value on number of stages to avoid exploring
                        the full graph. Default to 300
  -b BATCH_SIZE, --batch_size BATCH_SIZE
                        Maximum breath size. If it is reached, leads to hybrid
                        exploration. Default to 10 ** 10
  --alpha ALPHA         FWER Control value. Defaul to 10 ** (-8). A value of
                        zero will means that alpha is picked automatically.
  -v, --verbose         Wether to be verbose. If not specified (default), it
                        will not be verbose.
  --kmax KMAX           Maximum value on k to avoid exploring the full graph.
                        Default to 10 ** 8
  --dfs                 If this option is set, the code will be run in DFS
  --save_int            If this option is set, the current list of subgraphs
                        will be saved for future reruns.
  --restart             If this option is set, the code assumes that there is
                        a restart.

This program will save the list of significant connected subgraphs as well as
the associated p-values in the specified folder. Note that if runtime is too
slow, either increase the value of alpha, decrease sMax or increase the number
of cores. Caldera is quite memory intensive so make sure to keep track of
memory usage if running on a large file.
```
