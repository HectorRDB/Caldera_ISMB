# CALDERA

This is the repo where we develop the new code for our work on bacterial and meta GWAS work. It builds on [DBGWAS](https://gitlab.com/leoisl/dbgwas), described in [Jaillard et al](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007758).

## Installation

To use this, clone from github and install via pip

```{bash}
git clone git@github.com:HectorRDB/CALDERA.git
python3 -m pip install -e CALDERA/.
```

## Example

Then you can run the code. If the output of DBGWAS is stored in __/default_folder/out__, you first need to modify the files to create back the majority variant file.

```{bash}
toMajor.py -l /default_folder/out
```

Then you can run CALDERA using

```{bash}
caldera-script -l /default_folder/out
```

See `caldera-script --help`  for more options.
