# This program transforms the output from Caldera into formats 
# compatible with dbgwas step 3
import numpy as np
import pandas as pd
import os
from scipy.stats import chi2

def remove_redundant_subgraphs(sigs, sigs_pvals, verbose):
    """
    remove all subgraphs that are fully inside another subgraph
    and not significant
    
    Arguments:
        sigs: the significant graphs
        sigs_pvals: the p-value of those graphs
        verbose: whether to be verbose
    """
    not_all_checked = True
    i = 0
    while not_all_checked:
        S = sigs[i]
        inside_another = [all(node in sigs[j].nodes for node in S.nodes) and 
                            sigs_pvals[i] <= sigs_pvals[j] for j in range(i + 1, len(sigs))]
        if any(inside_another):
            sigs = np.delete(sigs, i)
            sigs_pvals = np.delete(sigs_pvals, i)
        else:
            i = i +1
        not_all_checked = (i != len(sigs))
        # 
    if verbose:
        print("We found {n} significant subgraphs".format(n = len(sigs)),
              flush = True)
    sigs = sigs[np.argsort(sigs_pvals)]
    sigs_pvals = np.sort(sigs_pvals)
    return sigs, sigs_pvals


def create_files(loc, out, sigs, sigs_pvals, G, verbose):
    """
    Create all the files necessary from the step1 folder from DBGWAS and
    the output from caldera's script
    
    Arguments:
        loc: Where the step1 folder is located
        out: The output folder is located. 
        sigs: the significant graphs
        sigs_pvals: the p-value of those graphs
        G: the structure
        verbose: whether to be verbose or not.
    """
    if not os.path.isdir(out + "/step1/"):
        os.mkdir(out + "/step1/")
    if not os.path.isdir(out + "/step2/"):
        os.mkdir(out + "/step2/")
    # Remove redundant subgraphs
    sigs, sigs_pvals = remove_redundant_subgraphs(sigs, sigs_pvals, verbose)
    # Create the pattern file
    patterns = pd.DataFrame()
    patterns['ID'] = range(len(sigs))
    patterns['pvals'] = [chi2.sf(th, 1) for th in sigs_pvals]
    patterns['statistics'] = sigs_pvals
    patterns['Pheno0'] = [S.ys[G.Pheno == 0].sum() for S in sigs]
    patterns['Pheno1'] = [S.ys[G.Pheno == 1].sum() for S in sigs]
    patterns['PhenoNA'] = [S.ys[np.isnan(G.Pheno)].sum() for S in sigs]
    patterns.to_csv(out + "/step2/patterns.txt", header = True,
                    index = None, sep = ' ', mode='a')
    # Create the gemma file
    # Get all the nodes that are in any of the significant subgraphs
    nodes = np.unique(np.concatenate([list(S.nodes) for S in sigs]))
    node_ccs = []
    for node in nodes:
        ccs_contain_nodes = [i for i in range(len(sigs)) if node in sigs[i].nodes]
        ccs_id = max(ccs_contain_nodes)
        node_ccs.append([node, ccs_id])
    nodes = np.array(node_ccs)
    others = np.array([node for node in range(G.pattern.shape[0]) if node not in nodes[:,0]])
    others = np.array([others, others + nodes.shape[0]], dtype = np.int64).transpose()
    all_nodes = np.concatenate([nodes, others])
    np.savetxt(out + "/step1/gemma_input.unitig_to_pattern.binary", all_nodes, fmt='%i')
    # Create the bugwas file
    with open(out + "/step1/bugwas_input.unique_rows_to_all_rows.binary", 'w') as f:
        for S in sigs:
            s = " ".join(map(str, list(S.nodes)))
            _ = f.write(s + "\n")
    # Create the transpose
    with open(out + "/step2/nodes_to_css", 'w') as f:
        for node in range(G.pattern.shape[0]):
            ccs_contain_nodes = [i for i in range(len(sigs)) if node in sigs[i].nodes]
            if (len(ccs_contain_nodes) == 0):
                ccs_contain_nodes = [-1]
            s = " ".join(map(str, ccs_contain_nodes))
            _ = f.write(s + "\n")
    # Create the other files
    weight = np.ones((G.pattern.shape[0], 1), dtype = np.int)
    np.savetxt(out + "/step1/weight_correction", weight, fmt='%i')
    os.symlink(loc + "/step1/graph.nodes", out + "/step1/graph.nodes")
    os.symlink(loc + "/step1/unitigs2PhenoCounter", out + "/step1/unitigs2PhenoCounter")
    os.symlink(loc + "/step1/phenoCounter", out + "/step1/phenoCounter")
    os.symlink(loc + "/step1/graph.edges.dbg", out + "/step1/graph.edges.dbg")
    os.system("touch " + out + "/step2/bugwas_out_barplot_BayesianWald_PCs.png")
    os.system("touch " + out + "/step2/bugwas_out_SNPs_PC_manhattan.png")
    os.system("touch " + out + "/step2/bugwas_out_tree_branchescolouredbyPC.png")
    print("Saving the output to " + out, flush = True)
