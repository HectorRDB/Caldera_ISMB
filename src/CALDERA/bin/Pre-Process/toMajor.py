#!/usr/bin/env python3

##### Script to get the right U matrix as input #####
import sys, getopt
import numpy as np

def main(argv):
    loc = ''
    try:
        opts, args = getopt.getopt(argv,"hl:",["help", "location="])
    except getopt.GetoptError:
        print('toMajor.py -l <folderLocation>')
        sys.exit(2)
    if opts == []:
        print('toMajor.py -l <folderLocation>')
        sys.exit()
    for opt, arg in opts:
        if (opt == '-h') | (opt == '--help') | (opt == ''):
            print ('toMajor.py -l <folderLocation>')
            sys.exit()
        else:
            loc = arg
    U = np.loadtxt(loc + "/step1/bugwas_input.all_rows.binary", skiprows = 1)
    U = U[:, 1:]
    flip = np.loadtxt(loc + "/step1/weight_correction")
    flip.shape=(-1, 1)
    # Re-transform everything back to original allele instead of minority allele
    U = (((U * 2 - 1) * flip + 1) / 2).astype(int)
    np.savetxt(loc + "/step1/major_pattern.all_rows.binary", U, fmt = "%i",
               delimiter = " ")

if __name__ == "__main__":
    main(sys.argv[1:])
