#!/usr/bin/env python3

import os, argparse, sys
import numpy as np
sys.path.append(os.path.abspath("/scr/AllUnitigs/"))
from OneStage import test_all_unitigs

parser = argparse.ArgumentParser(description =  "Test all unitigs")
parser.add_argument("-l", "--loc", help = "where to find the step 1 folder from DGWAS.",
                    type = str, required = True)
parser.add_argument("-o", "--output",
                    help = "Where to store the output. Default to loc/step2/",
                    type = str, default = "")
parser.add_argument("-C", "--communities",
                    help = "Number of communities to find. Default to 3",
                    type = int, default = 3)
parser.add_argument("-P", "--comFile",
                    help = "Location of the community assignement file." +
                            " If none is specified, default to k-means " +
                            "with k = communities",
                    type = str, default = "None")
parser.add_argument("--alpha",
                    help = "FWER Control value. Defaul to 10 ** (-8). A value of " +
                            "zero will means that alpha is picked automatically.",
                    type = np.float64, default = 10 ** (-8))
parser.add_argument("-v", "--verbose", dest = 'verbose', action = 'store_true',
                    help = "Wether to be verbose. If not specified (default), " +
                           "it will not be verbose.")
parser.set_defaults(verbose = False)
parser.parse_args()

args = parser.parse_args()
if args.output == "":
    output = args.loc + '/step2/'
else:
    output = args.output
if not os.path.exists(output):
    os.makedirs(output)
if args.comFile != "None":
    comFile = args.comFile
else:
    comFile = "None"
if args.verbose and comFile != "None":
    print("Reading community file from {file}".format(file = comFile), flush = True)
    print(flush = True)

alpha = args.alpha
if alpha <= 0:
    alpha = None
# Testing all unitigs

test_all_unitigs(loc = args.loc, output = output, communities = args.communities,
                 comFile = comFile, alpha = alpha, verbose = args.verbose)

# import os
# loc='/Data/Amikacin'
# alpha_str='0.000000'
# for power in range(1, 3):
#     alpha = 10 ** (-power)
#     print(alpha_str + "1")
#     out='Output/Amikacin/Caldera_S7_' + alpha_str + "1"
#     if not os.path.isdir(out + "/step1/"):
#         os.mkdir(out + "/step1/")
#     if not os.path.isdir(out + "/step2/"):
#         os.mkdir(out + "/step2/")
#     os.symlink(loc + "/step1/graph.nodes", out + "/step1/graph.nodes")
#     os.symlink(loc + "/step1/unitigs2PhenoCounter", out + "/step1/unitigs2PhenoCounter")
#     os.symlink(loc + "/step1/phenoCounter", out + "/step1/phenoCounter")
#     os.symlink(loc + "/step1/graph.edges.dbg", out + "/step1/graph.edges.dbg")
#     os.system("touch " + out + "/step2/bugwas_out_barplot_BayesianWald_PCs.png")
#     os.system("touch " + out + "/step2/bugwas_out_SNPs_PC_manhattan.png")
#     os.system("touch " + out + "/step2/bugwas_out_tree_branchescolouredbyPC.png")
#     alpha_str = alpha_str + '0'


# import os
# loc='/Data/Amikacin'
# alpha_str='0.00'
# for power in range(1, 8):
#     print(alpha_str + "1")
#     out='/Output/Amikacin/Caldera_S7_' + alpha_str + "1"
#     with open(out + '/step1/bugwas_input.unique_rows_to_all_rows.binary') as f:
#         lines = f.read().splitlines()
#     lines = [line.split(' ') for line in lines]
#     lines = [[int(i) for i in line] for line in lines]
#     with open(out + "/step2/nodes_to_css", 'w') as f:
#         for node in range(2356052):
#             ccs_contain_nodes = [i for i in range(len(lines)) if node in lines[i]]
#             if (len(ccs_contain_nodes) == 0):
#                 ccs_contain_nodes = [-1]
#             s = " ".join(map(str, ccs_contain_nodes))
#             _ = f.write(s + "\n")
#     alpha_str = alpha_str + '0'
