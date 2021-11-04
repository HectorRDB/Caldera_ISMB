#!/root/.pyenv/shims/python
################################################################################
################################################################################

import os, argparse
from random import normalvariate
import numpy as np
from Speed_Helper import *

parser = argparse.ArgumentParser(description =  "Test the speed of caldera vs coin script")
parser.add_argument("-i", "--id", help = "where to find the step 1 folder from DGWAS.",
                    type = int)
parser.add_argument("-t", "--timeout",
                    help = "Maximum time to run the script",
                    type = int, default = 2 * 24 * 60 * 60)
parser.add_argument("-n", "--number",
                    help = "Number of samples",
                    type = int, default = 100)
parser.add_argument("-C", "--Clusters",
                    help = "Number of populations",
                    type = int, default = 1)
parser.add_argument("--alpha",
                    help = "FWER Control value. Defaul to 0.05",
                    type = np.float64, default = 0.05)
parser.add_argument("-p", "--prop",
                    help = "Proportion of zeros and ones",
                    type = np.float64, default = 0.5)
parser.add_argument('--long', dest = 'long', action = 'store_true',
                    help = "Run the long simulation.")
parser.set_defaults(long = False)
parser.parse_args()
args = parser.parse_args()
# Run the code
random.seed(970)
np.random.seed(12)
if args.long:
    run_long_simu(id = args.id, timeout = args.timeout, n = args.number,
                  alpha = args.alpha, p = args.prop, C = args.Clusters)
else:
    run_simu(id = args.id, timeout = args.timeout, n = args.number,
            alpha = args.alpha, p = args.prop, C = args.Clusters)
