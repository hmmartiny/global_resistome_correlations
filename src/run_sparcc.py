
"""
SparCC code for GPU 
Adopted from https://github.com/bio-developer/sparcc/
"""

import os
#import pandas as pd
import cupy as cp
import argparse
import numpy as np

#from zero_replacement import zero_replacement
#from variation import variation_matrix


def parse_args():
    parser = argparse.ArgumentParser(
        "Run SparCC on GPU",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-d', '--data',
        type=str,
        required=True,
        help='File to data file with counts with parts as columns and samples as rows',
        dest='datafile'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output prefix or directory to store result files in',
        dest='output'
    )
    parser.add_argument(
        '--transpose',
        action='store_true',
        help="If true, transpose data to make parts as columns and samples as rows",
        dest='transpose'
    )
    parser.add_argument(
        '-n_iter', '--n_iterations',
        type=int,
        help='Number of SparCC iterations',
        default=50,
        dest='n_iter'
    )
    parser.add_argument(
        '-x_iter', '--n_exclusion_iterations',
        type=int,
        default=10,
        help='Number of exclusion iterations',
        dest='xiter'
    )
    parser.add_argument(
        '-n_perm', '--n_permutations',
        type=int,
        help='Number of bootstrapping samples to take to estimate p-values',
        default=5,
        dest='n_perm'
    )
    parser.add_argument(
        '-n_piter', '--n_perm_iterations',
        type=int,
        help='Number of iterations per permutation',
        default=5,
        dest='piter'
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Be verbose?',
        dest='verbose'
    )
    parser.add_argument(
        '-minFrags',
        type=int,
        default=50,
        help='Minimum frag value',
        dest='minFrags'
    )
    parser.add_argument(
        '-minSample',
        type=int,
        default=10,
        help='Minimum number of smaples with minFrag',
        dest='minSample'
    )
    parser.add_argument(
        '-tol',
        type=float,
        default=1e-3,
        dest='tol'
    )
    parser.add_argument(
        '-th',
        type=float,
        default=0.01,
        dest='th'
    )   

    return parser.parse_args()


def zero_replacement(matrix, alpha=0.5, n_samples=5000):
    counts = cp.zeros_like(matrix).T
    for i in range(0, matrix.shape[0]):
        p_matrix = cp.random.dirichlet(matrix[i, :] + alpha, size=(n_samples))
        counts[:, i] = p_matrix.mean(axis=0)
    return counts.T

def variation_matrix(matrix):
    """Calculate variations of matrix """
    matrix = cp.asarray(matrix)
    var_mat = cp.var(cp.log(matrix[:, :, None] * 1. / matrix[:, None]), axis=0)
    return var_mat
def variation_mat_slow(frame, shrink=False):
    '''
    Return the variation matrix of frame.
    Element i,j is the variance of the log ratio of components i and j.
    Slower version to be used in case the fast version runs out of memeory.
    '''
    frame_a = 1.* cp.asarray(frame)
    k    = frame_a.shape[1]
    V      = cp.zeros((k,k))
    for i in range(k-1):
        for j in range(i+1,k):
            y     = cp.array(cp.log(frame_a[:,i]/frame_a[:,j]))
            v = cp.var(y, ddof=1) # set ddof to divide by (n-1), rather than n, thus getting an unbiased estimator (rather than the ML one). 
            V[i,j] = v
            V[j,i] = v
    return V

def variation_mat(frame):
    '''
    Return the variation matrix of frame.
    Element i,j is the variance of the log ratio of components i and j.
    '''
    x    = 1.*cp.asarray(frame)
    n,m  = x.shape
    if m > 1000:
        return variation_mat_slow(frame)
    else:
        xx   = cp.tile(x.reshape((n,m,1)) ,(1,1,m))
        xx_t = xx.transpose(0,2,1)
        try:
            l    = cp.log(1.*xx/xx_t)
            V    = l.var(axis=0, ddof=1)
            return V
        except MemoryError:
            return variation_mat_slow(frame)
def basis_var(f, Var_mat, M, **kwargs):

    # compute basis variance
    M_inv = cp.linalg.inv(M)

    # sum elemenets along axis 1 to get t_i's
    V_vec = Var_mat.sum(axis=1)
    V_base = cp.dot(M_inv, V_vec)

    # if any are variances are less than zero, set them to V_min
    V_min = kwargs.get('V_min', 1e-10)
    V_base[V_base <= 0] = V_min

    return V_base

def C_from_V(Var_mat, V_base):
    Vi, Vj = cp.meshgrid(V_base, V_base)
    Cov_base = 0.5 * (Vi + Vj - Var_mat)
    C_base = Cov_base / cp.sqrt(Vi) / cp.sqrt(Vj)
    return C_base, Cov_base

def exclude_pair(C, th=0.1, previously_excluded=[]):
    C_temp = cp.triu(cp.abs(C), 1) # work only on upper triangle
    
    if len(previously_excluded) > 0:
        C_temp[tuple(cp.array(previously_excluded).T)] = 0
    i, j = cp.unravel_index(cp.argmax(C_temp), C_temp.shape)
    cmax = C_temp[i, j]
    if cmax > th:
        return i, j
    else:
        return None    

def run_sparcc(f, **kwargs):

    # setup parameters
    th = kwargs.get('th', 0.1) # exclusion threshold
    xiter = kwargs.get('xiter', 10) # exlusion iterations

    ## calculate variation matrix
    var_mat = variation_mat_slow(f)#variation_mat(f)
    var_mat_temp = var_mat.copy()
    
    # number of components
    D = var_mat.shape[0]
    M = cp.ones(shape=(D,D)) + cp.diag([D-2]*D)

    # approximate basis variance
    V_base = basis_var(f, var_mat_temp, M)
    C_base, Cov_base = C_from_V(var_mat, V_base)

    # iterate to exclude strongly correlated pairs
    excluded_pairs, excluded_comp = [], cp.array([])
    for xi in range(xiter):
        # search for a new pair to exclude
        to_exclude = exclude_pair(C_base, th, excluded_pairs)

        if to_exclude is None:
            break
        
        excluded_pairs.append(to_exclude)
        i, j = to_exclude
        M[i, j] -= 1
        M[j, i] -= 1
        M[i, i] -= 1
        M[j, j] -= 1
        inds = tuple(cp.array(excluded_pairs).T)
        var_mat_temp[inds] = 0
        var_mat_temp.T[inds] = 0

        # search for new
        nexcluded = cp.bincount(cp.ravel(cp.array(excluded_pairs)))
        excluded_comp_prev = set(excluded_comp.copy())
        excluded_comp = cp.where(nexcluded >= (D-3))[0]
        excluded_comp_new = set(excluded_comp) - excluded_comp_prev

        if len(excluded_comp_new) > 0:
            print(excluded_comp)
            if len(excluded_comp) > (D-4):
                print("Too many components excluded..")
                return
            
            for xcomp in excluded_comp_new:
                var_mat_temp[xcomp, :] = 0
                var_mat_temp[:, xcomp] = 0
                M[xcomp, :] = 0
                M[:, xcomp] = 0
                M[xcomp, xcomp] = 1
        
        # run another
        V_base = basis_var(f, var_mat_temp, M)
        C_base, Cov_base = C_from_V(var_mat, V_base)

        # set excluded components inferred values to nans
        for xcomp in excluded_comp:
            V_base[xcomp] = cp.nan
            C_base[xcomp, :] = cp.nan
            C_base[:, xcomp] = cp.nan
            Cov_base[xcomp, :] = cp.nan
            Cov_base[:, xcomp] = cp.nan
            
    return V_base, C_base, Cov_base

def basis_corr(f, **kwargs):

    k = f.shape[1] # number of components
    assert k >= 4, f'Not possible to detect correlation between composition of <4 components ({k} given)'

    tol = kwargs.get('tol', 1e-3) #tol = 1e-3 # tolerance for correlation range
    # actually run sparcc
    V_base, C_base, Cov_base = run_sparcc(f, **kwargs)
    if cp.max(cp.abs(C_base)) > 1 + tol:
        raise AssertionError("Sparsity assumption violated ({} > {}).".format(cp.max(cp.abs(C_base)), 1 + tol))
    return V_base, C_base, Cov_base

def sparcc(matrix, **kwargs):

    # put on gpu
    X = cp.asarray(matrix)

    # setup parameters
    n_iter = kwargs.get('iter', 5) # iterations to run
    verbose = kwargs.get('verbose', True) # print statements
    
    # list of outputs from different iterations 
    cor_list, var_list = [], []

    for i in range(n_iter):
        if verbose:
            print(f" Iteration: {i+1} / {n_iter}    ", end='\r')

        # convert to fractions with bayesian zero replacement
        frac = zero_replacement(X)

        # estimate basis correlations
        v_sparse, cor_sparse, cov_sparse = basis_corr(frac, **kwargs)
        var_list.append(cp.diag(cov_sparse))
        cor_list.append(cor_sparse)
    
    cor_array = cp.array(cor_list)
    var_array = cp.array(var_list)
    var_med = cp.nanmedian(var_array, axis=0) # median variances
    cor_med = cp.nanmedian(cor_array, axis=0) # median correlations
    x, y = cp.meshgrid(var_med, var_med)
    cov_med = cor_med * x**0.5 * y**0.
    
    return cor_med, cov_med

def get_axis(axis_name):

    if axis_name in [0, '0', 'rows', 'index']:
        return 0
    elif axis_name in [1, '1', 'cols', 'columns']:
        return 1
    else:
        return None

def permute_w_replacement(frame, axis=0):

    frame = cp.asarray(frame)

    
    s = frame.shape[axis]
    fun = lambda x: x[cp.random.randint(0, s, (1, s))][0]
    perm = cp.apply_along_axis(
        func1d = fun,
        axis=axis,
        arr = frame
    )
    return perm

def make_bootstraps(counts, nperm):

    counts = cp.asarray(counts)
    
    perms = {}
    for i in range(nperm):
        counts_perm = permute_w_replacement(counts)
    
        perms[i] = counts_perm
    
    return perms

def compare2sided(perm, real):
    return cp.abs(perm) >= cp.abs(real)

def compare1sided(perm, real):
    inds_abs = compare2sided(perm, real)
    inds_sign = cp.sign(perm) == cp.sign(real)
    return inds_abs & inds_sign

def get_pvalues(true_corrs, perm_corrs, test_type='two_sided', **kwargs):

    cmpfuns = {'two_sided': compare2sided, 'one_sided': compare1sided}
    cmpfun = cmpfuns[test_type]

    n_sig = cp.zeros(shape=true_corrs.shape)
    for perm_corr in perm_corrs:
        n_sig[cmpfun(perm_corr, true_corrs)] += 1
    
    n_perm = len(perm_corrs)
    p_vals = 1. * n_sig / n_perm
    p_vals[cp.diag_indices_from(p_vals)] = 1
    return p_vals

def run_bootstrapping(counts, **kwargs):

    nperm = kwargs.get('n_perm', 5)
    verbose = kwargs.get('verbose', True)

    # get the "real" correlations
    if verbose:
        print("Iterate on unshuffled data:")
    real_cor, real_cov = sparcc(counts, **kwargs)
    if verbose:
        print("Finished iteration on unshuffled data.")

    bootstrap_corrs = []
    kwargs['iter'] = kwargs['piter']
    for i in range(nperm):
        if verbose:
            print(f"Bootstrap iteration: {i+1} / {nperm}   ", end="\r")
        # make permutation
        counts_perm = permute_w_replacement(counts)

        # run sparcc on shuffled data
        perm_cor, perm_cov = sparcc(counts_perm, **kwargs)
        bootstrap_corrs.append(perm_cor)

    
    test_type = kwargs.get('test_type', 'two_sided')
    p_values = get_pvalues(real_cor, bootstrap_corrs, test_type) 
    return p_values, real_cor, real_cov

if __name__ == "__main__":
    args = parse_args()
    if args.verbose:
        s = "Parameter set as follows\n"
        s += f"* Number of iterations: {args.n_iter}\n"
        s += f"* Number of exclusion iterations: {args.xiter}\n"
        s += f"* Number of permutations: {args.n_perm}\n"
        s += f"* Number of permutation iterations: {args.piter}\n"
        print(s)
    
    # count number of columns first
    with open(args.datafile, 'r') as f: firstLine = f.readline()
    cols = firstLine.strip().split(',')
    nCols = len(cols)

    df = cp.loadtxt(args.datafile, usecols=tuple(range(1, nCols)), skiprows=1, delimiter=',')
#    df = pd.read_csv(args.datafile, sep='\t' if args.datafile.endswith('.tsv') else ',', index_col=0)
    if args.verbose:
        print("Loaded data from file:", args.datafile)
        print(" shape:", df.shape)

    # drop columns that either only contain zeroes or 
    # where there are too few samples with enough fragments 
    minFrags = args.minFrags
    minSample = args.minSample  #df.shape[0] * minSamplePCT
    if args.verbose:
        print("Applying following filters:")
        print(" * Minimum fragment count per sample for a gene:", minFrags)
        print(" * at least {} should have minFrags".format(minSample)) 
    
    #nonzero_idx = cp.where(cp.any(df > 0, axis=0))
    nonzero_idx = cp.where(cp.sum(df >= minFrags, axis=0) >= minSample) #cp.where(cp.sum(df >= minFrags, axis=1) >= minSample)
    df = df[:, nonzero_idx[0]]
    if args.verbose:
        print("After filtering, data shape is:", df.shape)
    if df.shape[1] == 0:
        raise ValueError("No column passed filters, exiting..")

    # save columns that did not contain zeros
    dataCols = np.asarray(cols[1:])
    nonzeroCols = dataCols[nonzero_idx[0].get()]
    np.save(file=args.output + 'columns', arr=nonzeroCols)

#    if args.transpose:
#        df = df.reset_index().drop_duplicates(subset=['run_accession', 'refSequence']).pivot(
#            index='run_accession', columns='refSequence', values='fragmentCountAln'
#        ).fillna(0)
#        if args.verbose: print("Transposed data")

    # start bootstrapping
    print("Running SparCC bootstrapping procedure..")
    p_values, real_corr, real_cov = run_bootstrapping(
        df, 
        n_perm=args.n_perm,
        iter=args.n_iter,
        xiter=args.xiter,
        piter=args.piter,
        verbose=args.verbose,
        tol=args.tol,
        th=args.th
    )
    print("Running SparCC bootstrapping procedure.. done")

    # save
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    if args.verbose:
        print("Writing to:", args.output)
    
    # save arrays as binary files
    cp.save(file=args.output + 'sparr_corr', arr=real_corr)
    cp.save(file=args.output + 'sparr_pval', arr=p_values)


#    # write correlations
#    corrDF = pd.DataFrame(real_corr.get(), columns=df.columns, index=df.columns)
#
#    # write pval
#    pvalDF = pd.DataFrame(p_values.get(), columns=df.columns, index=df.columns)
#    try:
#        corrDF.to_csv(args.output + 'sparcc_corr.csv', float_format="%.5f")
#        pvalDF.to_csv(args.output + 'sparcc_pval.csv', float_format="%.15f")
#    except TypeError:
#        print("Failed to save to csv")
