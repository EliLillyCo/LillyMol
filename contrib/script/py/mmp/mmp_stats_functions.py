###################################################################
""" Summary: Useful stats functions for working with MMP data

"""
###################################################################
import numpy as np
import pandas as pd

from math import sqrt, log10
from scipy.stats import t

# things needed for unit testing only
import unittest
import tempfile


############################################################ 
#
# Custom aggregation functions
# 
# See unit tests for usage as well as mmp_pairs_objects for usage
# 
############################################################

def diffn(x, level, inc_low_n_vals=False):
    """This function was developed for use in
    aggregating certain MMP data such as Microsomal Metabolic Turnover (diff60)
    or Single Point CYPxxx Inhibition (Diff50 of Est pIC50). Note that it will
    return Nan for array with <3 values"""

    n_ = len(x)

    # timing data inc_low_n_vals=False
    # ncalls  tottime  percall  cumtime  percall
    # 20000    0.169    0.000    7.672    0.000
    # 20000    0.170    0.000    7.569    0.000
    # 20000    0.171    0.000    7.733    0.000
    # 20000    0.169    0.000    7.752    0.000
    #
    #
    if n_ < 3 and inc_low_n_vals is False:
        return np.nan, np.nan, np.nan

    # Baseline with full diffn calc on n<3
    # ncalls  tottime  percall  cumtime  percall
    # 20000    0.339    0.000   13.810    0.001
    # 20000    0.344    0.000   13.591    0.001
    # 20000    0.333    0.000   13.691    0.001
    # 20000    0.320    0.000   13.551    0.001
    #
    # adding below catch for n==1 and inc_low_n_vals is True
    # ncalls  tottime  percall  cumtime  percall
    # 20000    0.257    0.000   10.829    0.001
    # 20000    0.257    0.000   10.948    0.001
    # 20000    0.268    0.000   10.984    0.001
    # 20000    0.261    0.000   11.072    0.001
    #
    # ----------------------------------------
    # 2019-04-05 Baseline (hardware changes)
    # ncalls  tottime  percall  cumtime  percall
    # 20000    0.165    0.000    6.151    0.000
    # 20000    0.159    0.000    5.890    0.000
    # 20000    0.159    0.000    6.053    0.000
    # 20000    0.174    0.000    6.673    0.000
    #
    # edits to add round(x, 4) to all returns
    #
    # ncalls  tottime  percall  cumtime  percall
    # 20000    0.176    0.000    6.191    0.000
    # 20000    0.179    0.000    6.470    0.000
    # 20000    0.176    0.000    6.100    0.000
    # 20000    0.179    0.000    6.243    0.000
    #
    elif n_ == 1 and inc_low_n_vals is True:
        #
        avg_ = np.average(x)
        lx_ = log10(level / (100.0 - level))
        a_ = (10.0 ** (avg_ + lx_))
        diffn = round(100.0 * a_ / (1.0 + a_) - level, 4)

        return np.nan, diffn, np.nan

    avg_ = np.average(x)
    sd_ = np.std(x, ddof=1)
    interval_ = t.interval(0.9, df=(n_ - 1))[1] * sd_ / sqrt(n_)
    low90_ = avg_ - interval_
    upp90_ = avg_ + interval_
    lx_ = log10(level / (100.0 - level))

    try:
        a_ = (10.0 ** (avg_ + lx_))
        diffn = round(100.0 * a_ / (1.0 + a_) - level, 4)
    except:
        diffn = np.nan

    try:
        b_ = (10.0 ** (low90_ + lx_))
        diffn_low = round(100.0 * b_ / (1.0 + b_) - level, 4)
    except:
        diffn_low = np.nan

    try:
        c_ = (10.0 ** (upp90_ + lx_))
        diffn_upp = round(100.0 * c_ / (1.0 + c_) - level, 4)
    except:
        diffn_upp = np.nan

    # convert to pandas Series to allow return of column labels (index)
    # worked once! but not twice :-)
    # s_ = pd.Series([diffn, diffn_upp, diffn_low], index=['diffn', 'diffn_upp', 'diffn_low'])

    return diffn_low, diffn, diffn_upp


def diffn_list_rtn(n, inc_low_n_vals):
    """returns a list containing: 
    [diffn, diffn_upp, diffn_low]"""
    def diff_(x):
        # note below returns full list: [diffn, diffn_upp, diffn_low]
        return diffn(x, n, inc_low_n_vals)
    diff_.__name__ = 'diff%s_all' % n
    return diff_


def diffn_agg(n):
    """returns only diffn"""
    def diff_(x):
        return diffn(x, n)[1]
    diff_.__name__ = 'diff%s' % n
    return diff_  


def diffn_agg_upp(n):
    """Returns only diffn_upp"""
    def diff_(x):
        return diffn(x, n)[2]
    diff_.__name__ = 'diff%s_upp' % n
    return diff_  


def diffn_agg_low(n):
    """Returns only diffn_upp"""
    def diff_(x):
        return diffn(x, n)[0]
    diff_.__name__ = 'diff%s_low' % n
    return diff_  


#####################################################
#
# An alternate aggregation function for use with data such as 
# Unbound Microsomal Intrinsic Clearance (rClint-u) where
# we want to calculate a 'Fold Change' which is the mean diff 
# in log scale for Pair_L/Pair_R => Index.  Additionally the 
# result becomes -10 if negative 0.1 change)
#
#####################################################

def _mean_diff_invlog(x, inc_low_n_vals=False):
    """Aggregate function for use with deltas in log scale such as
    Unbound Microsomal Intrinsic Clearance (rClint-u) where we want
    to calculate a 'Fold Change' which is the mean diff in log
    scale for Pair_L/Pair_R => Index.  Additionally the result
    becomes -10 if negative 0.1 change.  It is expected that the
    input is log transformed but the output will be transformed back
    into original units via inverse log transform.  Equivalent to the
    Mean ratio or Fold Change or Mean difference
    """
    
    n_ = len(x)

    # These two lines have been speed tested with external code found in mmp_stats_functions_timer.py
    # This function was run 20,000 times on lists of randomly generated but fixed seed numbers. Of these
    # 20K lists, 10K are of length <3 and rest are of length between 3 to 100. Execution time is measured
    # by cProfile and the cumtime (cumulative total) is reported below. Run 4 times to smooth out cpu peaks.
    # The data is now identical in each run so we can compare times effectively/reproducibly.
    #
    # (1) Without 'CRITICAL TWO LINES':
    # ncalls	tottime	percall	cumtime	percall
    # 20000	0.269	0	13.744	0.001
    # 20000	0.262	0	13.668	0.001
    # 20000	0.265	0	13.546	0.001
    # 20000	0.273	0	13.662	0.001
    #
    # (2) With line "if n_ < 3: return np.nan, np.nan, np.nan":
    # ncalls	tottime	percall	cumtime	percall
    # 20000	0.139	0	7.604	0
    # 20000	0.137	0	7.598	0
    # 20000	0.137	0	7.549	0
    # 20000	0.14	0	7.599	0
    #
    # (3)(a) With line "if n_ < 3 and inc_low_n_vals is False: return np.nan, np.nan, np.nan":
    # and ensuring call with inc_low_n_vals=False
    # ncalls	tottime	percall	cumtime	percall
    # 20000    0.148    0.000    8.060    0.000
    # 20000    0.141    0.000    7.561    0.000
    # 20000    0.141    0.000    7.654    0.000
    # 20000    0.152    0.000    7.717    0.000
    #
    # (3)(b) With line "if n_ < 3 and inc_low_n_vals is False: return np.nan, np.nan, np.nan":
    # and ensuring call with inc_low_n_vals=True  <----- *NB*
    # ncalls	tottime	percall	cumtime	percall
    # 20000    0.270    0.000   13.642    0.001
    # 20000    0.273    0.000   13.872    0.001
    # 20000    0.275    0.000   13.919    0.001
    # 20000    0.266    0.000   13.572    0.001
    #
    # CRITICAL 2 LINES
    # if n_ < 3 and inc_low_n_vals is False:
    #    return np.nan, np.nan, np.nan

    # having optimised the above and chosen option (3) we now add this line to catch n_=1 cases
    # and test call with inc_low_n_vals=True --> The results are good wso we'll leave the line in
    # ncalls	tottime	percall	cumtime	percall
    # 20000    0.195    0.000   10.806    0.001
    # 20000    0.195    0.000   10.739    0.001
    # 20000    0.198    0.000   10.871    0.001
    # 20000    0.216    0.000   10.850    0.001
    #
    # and test call with inc_low_n_vals=False
    # 20000    0.142    0.000    7.757    0.000
    # 20000    0.143    0.000    7.756    0.000
    # 20000    0.150    0.000    8.116    0.000
    # 20000    0.145    0.000    7.787    0.000
    #
    # elif n_ == 1 and inc_low_n_vals is True:
    #    return np.nan, x[0], np.nan

    # and finally, having fixed the above two lines and confirmed that they do improve performance (speed)
    # test that the below change also helps
    # test call with inc_low_n_vals=True --> The results are good wso we'll leave the line in
    # ncalls	tottime	percall	cumtime	percall
    # 20000    0.156    0.000    7.839    0.000
    # 20000    0.150    0.000    7.858    0.000
    # 20000    0.150    0.000    7.689    0.000
    # 20000    0.148    0.000    7.864    0.000
    #
    # test call with inc_low_n_vals=False
    # ncalls	tottime	percall	cumtime	percall
    # 20000    0.147    0.000    7.756    0.000
    # 20000    0.143    0.000    7.573    0.000
    # 20000    0.142    0.000    7.638    0.000
    # 20000    0.143    0.000    7.621    0.000
    #
    # elif n_ == 2 and inc_low_n_vals is True:
    #    return np.nan, np.average(x), np.nan

    # Here's a version that actually works.  Previous one fails
    # as can't address x as array so x[0] fails
    #
    # ncalls	tottime	percall	cumtime	percall
    # 20000    0.151    0.000    7.798    0.000
    # 20000    0.158    0.000    8.155    0.000
    # 20000    0.153    0.000    7.887    0.000
    # 20000    0.149    0.000    7.726    0.000
    #
    # ncalls	tottime	percall	cumtime	percall
    # 20000    0.144    0.000    7.745    0.000
    # 20000    0.141    0.000    7.537    0.000
    # 20000    0.145    0.000    7.822    0.000
    # 20000    0.139    0.000    7.618    0.000
    #
    if n_ < 3:
        if inc_low_n_vals is False:
            return np.nan, np.nan, np.nan
        else:
            return np.nan, (10 ** np.average(x)), np.nan
    
    # Mean of these deltas in log scale
    avg_ = np.average(x)
    # sd of these deltas in log scale
    sd_ = np.std(x, ddof=1)
    int_ = t.interval(0.9, df=(n_ - 1))[1]

    low90_ = avg_ - int_ * (sd_ / sqrt(n_))
    upp90_ = avg_ + int_ * (sd_ / sqrt(n_))
    
    # need to reverse log:
    fold_change = round((10 ** avg_), 4)
    fold_change_upper90 = round((10 ** upp90_), 4)
    fold_change_lower90 = round((10 ** low90_), 4)
    
    return fold_change_lower90, fold_change, fold_change_upper90


def _mean_diff(x, inc_low_n_vals=False):
    """Same as fold_change but no inverse log function
    after aggregation.  For use with data like LogP/D
    """

    n_ = len(x)

    # see timing data for _mean_diff_inv_log version
    if n_ < 3:
        if inc_low_n_vals is False:
            return np.nan, np.nan, np.nan
        else:
            return np.nan, np.average(x), np.nan

    # Mean of these deltas in log scale
    avg_ = round(np.average(x), 4)
    # sd of these deltas in log scale
    sd_ = np.std(x, ddof=1)
    int_ = t.interval(0.9, df=(n_ - 1))[1]

    low90_ = round(avg_ - int_ * (sd_ / sqrt(n_)), 4)
    upp90_ = round(avg_ + int_ * (sd_ / sqrt(n_)), 4)

    return low90_, avg_, upp90_


# Added these to allow the previous mean_diff and mean_diff_invlog functions to be used
# with new keyword args, specifically to print out all low n value results. Returning the
# low n value results slows the function down and they are not useful in some MMP analysis
# but in other cases such as MMP prediction work we still need these data points.
def mean_diff(inc_low_n_vals):
    def mean_diff_(x):
        # note below returns full list: [diffn, diffn_upp, diffn_low]
        return _mean_diff(x, inc_low_n_vals)
    mean_diff_.__name__ = 'mean_diff'
    return mean_diff_


def mean_diff_invlog(inc_low_n_vals):
    def mean_diff_invlog_(x):
        # note below returns full list: [diffn, diffn_upp, diffn_low]
        return _mean_diff_invlog(x, inc_low_n_vals)
    mean_diff_invlog_.__name__ = 'mean_diff_invlog'
    return mean_diff_invlog_


#####################################################
#
# More useful stuff for counting number of positive
# or negative values in df
#
#####################################################

def n_neg_diff(x):
    return sum(n < 0 for n in x)


def n_pos_diff(x):
    return sum(n > 0 for n in x)


#####################################################
#
# Functions to use with categorical data
#
#####################################################
def category_moved_up_one(x):
    return sum(n == 1 for n in x)


def category_moved_up_one_pct(x):
    return round(float(sum(n == 1 for n in x))/float(len(x))*100, 4)
    
    
def category_moved_down_one(x):
    return sum(n == -1 for n in x)
    
    
def category_moved_down_one_pct(x):
    return round(float(sum(n == -1 for n in x))/float(len(x))*100, 4)


def category_no_change(x):
    return sum(n == 0 for n in x)


def category_no_change_pct(x):
    return round(float(sum(n == 0 for n in x))/float(len(x))*100, 4)


def category_changed_class(x):
    return sum(n != 0 for n in x)


def category_changed_class_pct(x):
    return round(float(sum(n != 0 for n in x))/float(len(x))*100, 4)


#####################################################
#
# unittest everything
#
#####################################################

class _TestMMPStatsFunctions(unittest.TestCase):
    """Test class to test the object and methods"""
    
    def setUp(self):
        #
        self.maxDiff = None
        
        self.test_data_diff = [-0.542, -3.043, 0.264, 0.094, -0.262, 0.344, 0.769, 0.811, -1.350, -1.475, -0.027,
                               -3.200, -3.345, -2.950]

        self.test_dataframe = pd.DataFrame({
            'A': ['aa', 'ab', 'ac', 'ad', 'aa', 'ab', 'ac', 'ad', 'aa', 'aa', 'ab'],
            'B': [-0.542, -3.043, 0.264, 0.094, -0.262, 0.344, 0.094, -0.262, -0.555, -0.54, -0.27],
            'C': [-3.043, 0.264, 0.094, -0.262, 0.344, 0.769, 0.094, -0.262, -3.001, -3.10, 0.35],
            'D': [0.264, 0.094, -0.262, 0.344, 0.769, 0.811, 0.094, -0.262, -0.260, -3.001, 0.100],
            'E': [-0.262, 0.344, 0.769, 0.811, -1.350, -1.475, 0.094, -0.262, -0.254, -0.254, 0.901],
            'F': [-0.262, 0.344, 0.769, 0.811, -1.350, -1.475, 0.094, -0.262, -0.256, -0.206, 0.344]
        })
        
        self.test_data_diff_goldresult_diffn = -46.791482
        self.test_data_diff_goldresult_diffn_upp = -14.72599757
        self.test_data_diff_goldresult_diffn_low = -57.27662627
       
        self.test_dataframe_goldenresult = pd.DataFrame({
            'A': ['aa', 'ab', 'ac', 'ad'],
            'C_count': [4, 3, 2, 2],
            'C_mean': [-2.200, 0.461, 0.094, -0.262],
            'C_diff60': [-59.0624, 21.2594, np.nan, np.nan],
            'C_diff60_low': [-59.9905, 0.3045, np.nan, np.nan],
            'C_diff60_upp': [-11.5932, 32.5238, np.nan, np.nan]
        })
        self.test_dataframe_goldenresult.sort_index(axis=1, inplace=True)
       
        self.test_diff_2_inputfile = tempfile.NamedTemporaryFile(delete=False)

        self.test_diff_2_input = pd.DataFrame({
            'FRAG_R': ['O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O=[1CH]C',
                       'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', '[1CH3]N1CCO[2CH2]C1'],
            'FRAG_L': ['s1[1cH][n]cc1', 's1[1cH][n]cc1', 's1[1cH][n]cc1', 's1[1cH][n]cc1', 's1[1cH][n]cc1',
                       's1[1cH][n]cc1', 's1[1cH][n]cc1', 's1[1cH][n]cc1', 's1[1cH][n]cc1', 's1[1cH][n]cc1',
                       's1[1cH][n]cc1', 's1[1cH][n]cc1', 'O=[1CH]N1C[2CH2]CCC1'],
            'logitM_DIFF': [-3.380, -2.895, -3.380, -1.282, -2.170, -0.226, -1.291, -1.715, -0.798, -1.378, -0.738,
                            -2.426, 7]
        })

        self.test_diff_2_goldenout = pd.DataFrame({
                'FRAG_L': ['s1[1cH][n]cc1', 's1[1cH][n]cc1', 'O=[1CH]N1C[2CH2]CCC1'],
                'FRAG_R': ['O[1CH]=O', 'O=[1CH]C', '[1CH3]N1CCO[2CH2]C1'],
                'diff60': [-59.108675, -54.253683, np.nan],
                'diff60_low90': [-59.918216, -58.174201, np.nan],
                'diff60_upper90': [-51.007039, -43.343025, np.nan],
                'num_vals': [6, 6, 1]
        })

        self.test_diff_2_goldenout_low_n_val = pd.DataFrame({
            'FRAG_L': ['s1[1cH][n]cc1', 's1[1cH][n]cc1', 'O=[1CH]N1C[2CH2]CCC1'],
            'FRAG_R': ['O[1CH]=O', 'O=[1CH]C', '[1CH3]N1CCO[2CH2]C1'],
            'diff60': [-59.108675, -54.253683, 7.00000],
            'diff60_low90': [-59.918216, -58.174201, np.nan],
            'diff60_upper90': [-51.007039, -43.343025, np.nan],
            'num_vals': [6, 6, 1]
        })

        self.test_foldchange_input = pd.DataFrame({
            'FRAG_L': ['O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O',
                       'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O',
                       'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O',
                       'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', '[1CH3]N1CCO[2CH2]C1', 'c1ccccc1', 'c1ccccc1'],
            'FRAG_R': ['O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C',
                       'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C',
                       'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C',
                       'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]N1C[2CH2]CCC1', 'c1ccccc1', 'c1ccccc1'],
            'pLog_rClintu_DIFF': [-1.99704061, -1.806422468, -1.827428534, -1.90152502, -1.249553199, -0.664484291,
                                  -1.857100559, -1.5845214, -2.667472973, -1.479999458, -1.635125098, -2.415917968,
                                  -1.823722494, -2.52674842, -1.418294469, -1.19056738, -1.955927983, -1.285681124,
                                  -1.27009994, -1.395897893, -1.641798956, -2.063255775, -1.529434468, -2.43015166,
                                  -1.620558471, -1.2, 0.2, 0.3]
        })

        # data from R: 77.4282849546501, 53.6474674395413,37.1705348292594
        self.test_foldchange_goldenout = pd.DataFrame({
            'FRAG_L': ['O[1CH]=O', '[1CH3]N1CCO[2CH2]C1', 'c1ccccc1'],
            'FRAG_R': ['O=[1CH]C', 'O=[1CH]N1C[2CH2]CCC1', 'c1ccccc1'],
            'pLog_rClintu_DIFF_count': [25, 1, 2],
            'pLog_rClintu_DIFF_mean_diff_invlog': [(0.0129, 0.0186, 0.0269),
                                                   (np.nan, np.nan, np.nan),
                                                   (np.nan, np.nan, np.nan)]
        })

        self.test_mean_diff_goldenout = pd.DataFrame({
            'FRAG_L': ['[1CH3]N1CCO[2CH2]C1'],
            'FRAG_R': ['O=[1CH]N1C[2CH2]CCC1'],
            'Solubility_pH2_mMolar_DIFF_count': [15],
            'Solubility_pH2_mMolar_DIFF_mean_diff': [(-2.4688, -1.9984, -1.528)]
        })

        self.test_foldchange_goldenout_inc_low_n_vals = pd.DataFrame({
            'FRAG_L': ['O[1CH]=O', '[1CH3]N1CCO[2CH2]C1'],
            'FRAG_R': ['O=[1CH]C', 'O=[1CH]N1C[2CH2]CCC1'],
            'pLog_rClintu_DIFF_count': [25, 1],
            'pLog_rClintu_DIFF_mean_diff_invlog': [(0.012915177, 0.018640209, 0.026903030),
                                                   (np.nan, -1.20000, np.nan)]
        })

        self.test_foldchange_input_2 = pd.DataFrame({
            'FRAG_L': ['[1CH3]N1CCO[2CH2]C1', '[1CH3]N1CCO[2CH2]C1', '[1CH3]N1CCO[2CH2]C1', '[1CH3]N1CCO[2CH2]C1',
                       '[1CH3]N1CCO[2CH2]C1', '[1CH3]N1CCO[2CH2]C1', '[1CH3]N1CCO[2CH2]C1', '[1CH3]N1CCO[2CH2]C1',
                       '[1CH3]N1CCO[2CH2]C1', '[1CH3]N1CCO[2CH2]C1', '[1CH3]N1CCO[2CH2]C1', '[1CH3]N1CCO[2CH2]C1',
                       '[1CH3]N1CCO[2CH2]C1', '[1CH3]N1CCO[2CH2]C1', '[1CH3]N1CCO[2CH2]C1'],
            'FRAG_R': ['O=[1CH]N1C[2CH2]CCC1', 'O=[1CH]N1C[2CH2]CCC1', 'O=[1CH]N1C[2CH2]CCC1', 'O=[1CH]N1C[2CH2]CCC1',
                       'O=[1CH]N1C[2CH2]CCC1', 'O=[1CH]N1C[2CH2]CCC1', 'O=[1CH]N1C[2CH2]CCC1', 'O=[1CH]N1C[2CH2]CCC1',
                       'O=[1CH]N1C[2CH2]CCC1', 'O=[1CH]N1C[2CH2]CCC1', 'O=[1CH]N1C[2CH2]CCC1', 'O=[1CH]N1C[2CH2]CCC1',
                       'O=[1CH]N1C[2CH2]CCC1', 'O=[1CH]N1C[2CH2]CCC1', 'O=[1CH]N1C[2CH2]CCC1'],
            'Solubility_pH2_mMolar_DIFF': [-0.02380751, -2.648036977, -0.855180236, -2.612320185, -2.732343491,
                                                -0.052453718, -2.240123595, -1.344902599, -2.532140934, -2.819401653,
                                                -2.545212373, -2.832473093, -2.920275198, -1.088678392, -2.728023456]
            # 'Solubility_pH2_mMolar_DIFF': [0.016217488, -2.030728258, -1.747948595, -1.598866354, -1.936057396, \
            #                                     0.045988318, -1.764970198, -1.948350161, -2.027001232, -2.029993245, \
            #                                     -2.020599475, -2.023591488, -2.013101327, -1.599220664, -2.070137256]
        })

        # data from R:
        self.test_foldchange_goldenout_2 = pd.DataFrame({
            'FRAG_L': ['[1CH3]N1CCO[2CH2]C1'],
            'FRAG_R': ['O=[1CH]N1C[2CH2]CCC1'],
            'Solubility_pH2_mMolar_DIFF_count': [15],
            'Solubility_pH2_mMolar_DIFF_mean_diff_invlog': [(0.0034, 0.01, 0.0297)]
        })

        # categorical data
        self.test_categorical_dataframe = pd.DataFrame({
            'FRAG_L': ['O[1CH]=O', 'O[1CH]=O', 'O[1CH]=O', '[1CH3]N1CCO[2CH2]C1', '[1CH3]N1CCO[2CH2]C1'],
            'FRAG_R': ['O=[1CH]C', 'O=[1CH]C', 'O=[1CH]C', 'O=[1CH]N1C[2CH2]CCC1', 'O=[1CH]N1C[2CH2]CCC1'],
            'CAT_DATA_DIFF': [0, 0, 1, 0, -1]
        })

        # categorical data golden
        self.test_categorical_dataframe_golden = pd.DataFrame({
            'FRAG_L': ['O[1CH]=O', '[1CH3]N1CCO[2CH2]C1'],
            'FRAG_R': ['O=[1CH]C', 'O=[1CH]N1C[2CH2]CCC1'],
            'CAT_DATA_DIFF_category_changed_class': [1, 1],
            'CAT_DATA_DIFF_category_moved_down_one': [0, 1],
            'CAT_DATA_DIFF_category_moved_up_one': [1, 0],
            'CAT_DATA_DIFF_category_no_change': [2, 1],
            'CAT_DATA_DIFF_count': [3, 2]
        })

        self.test_categorical_dataframe_pct_golden = pd.DataFrame({
            'CAT_DATA_DIFF_category_changed_class_pct': [33.3333, 50.0],
            'CAT_DATA_DIFF_category_moved_down_one_pct': [0, 50],
            'CAT_DATA_DIFF_category_moved_up_one_pct': [33.3333, 0],
            'CAT_DATA_DIFF_category_no_change_pct': [66.6667, 50.0],
            'CAT_DATA_DIFF_count': [3, 2],
            'FRAG_L': ['O[1CH]=O', '[1CH3]N1CCO[2CH2]C1'],
            'FRAG_R': ['O=[1CH]C', 'O=[1CH]N1C[2CH2]CCC1']
        })

    def test_diff_array(self):
        
        (result_diffn_low, result_diffn, result_diffn_upp) = diffn(self.test_data_diff, 60)

        # pp.pprint(result_diffn)
        # pp.pprint(result_diffn_upp)
        # pp.pprint(result_diffn_low)

        self.assertAlmostEqual(result_diffn, self.test_data_diff_goldresult_diffn, 4)
        self.assertAlmostEqual(result_diffn_upp, self.test_data_diff_goldresult_diffn_upp, 4)
        self.assertAlmostEqual(result_diffn_low, self.test_data_diff_goldresult_diffn_low, 4)

    def test_diff_df_individual(self):
        
        # group, affregate, convert object to df, sort index
        f = {'C': ['mean', 'count', diffn_agg_low(60), diffn_agg(60), diffn_agg_upp(60)]}
        grouped_stats = self.test_dataframe.groupby(['A']).agg(f)
        grouped_stats.columns = ['_'.join(col).strip() for col in grouped_stats.columns.values]
        # this is important for unittesting as order or columns can change and needs to be made consistent
        grouped_stats = pd.DataFrame(grouped_stats).reset_index()
        grouped_stats.sort_index(axis=1, inplace=True)
        
        # Use pandas assert on data frames
        # pp.pprint(self.test_dataframe_goldenresult)
        # pp.pprint(grouped_stats)
        pd.util.testing.assert_frame_equal(self.test_dataframe_goldenresult, grouped_stats, check_exact=False)

    def test_diff_df_series1(self):

        # group, affregate, convert object to df, sort index
        # print self.test_dataframe.head(10)
        f = {'C': ['mean', 'count', diffn_list_rtn(60, inc_low_n_vals=False)]}
        grouped_stats = self.test_dataframe.groupby(['A'], sort=False).agg(f)  # [diffn_list_rtn(60), len])
        grouped_stats.columns = ['_'.join(col).strip() for col in grouped_stats.columns.values]
        grouped_stats[['C_diff60_low', 'C_diff60', 'C_diff60_upp']] = grouped_stats['C_diff60_all'].apply(pd.Series)
        grouped_stats.drop('C_diff60_all', axis=1, inplace=True)
        # this is important for unittesting as order or columns can change and needs to be made consistent
        grouped_stats = pd.DataFrame(grouped_stats).reset_index()
        grouped_stats.sort_index(axis=1, inplace=True)
        # pp.pprint(self.test_dataframe_goldenresult)
        # pp.pprint(grouped_stats)
        pd.util.testing.assert_frame_equal(self.test_dataframe_goldenresult, grouped_stats)

    def test_diff_df_series2(self):

        f = {'logitM_DIFF': ['count', diffn_list_rtn(60, inc_low_n_vals=False)]}
        grouped_stats = self.test_diff_2_input.groupby(['FRAG_L', 'FRAG_R'], sort = False).agg(f)
        grouped_stats.columns = ['_'.join(col).strip() for col in grouped_stats.columns.values]
        grouped_stats[['diff60_low90', 'diff60', 'diff60_upper90']] = grouped_stats['logitM_DIFF_diff60_all'].apply(pd.Series)
        grouped_stats.drop('logitM_DIFF_diff60_all', axis=1, inplace=True)
        grouped_stats.rename(columns={'logitM_DIFF_count': 'num_vals'}, inplace=True)

        grouped_stats = pd.DataFrame(grouped_stats).reset_index()
        grouped_stats.sort_index(axis=1, inplace=True)

        pd.util.testing.assert_frame_equal(grouped_stats, self.test_diff_2_goldenout)

    def test_mean_diff_invlog_df(self):

        f = {'pLog_rClintu_DIFF': ['count', mean_diff_invlog(inc_low_n_vals=False)]}
        grouped_stats = self.test_foldchange_input.groupby(['FRAG_L', 'FRAG_R'], sort=False).agg(f)
        grouped_stats.columns = ['_'.join(col).strip() for col in grouped_stats.columns.values]
        
        grouped_stats = pd.DataFrame(grouped_stats).reset_index()
        grouped_stats.sort_index(axis=1, inplace=True)

        pd.util.testing.assert_frame_equal(grouped_stats, self.test_foldchange_goldenout)

    def test_mean_diff_invlog_df_2(self):
    
        f = {'Solubility_pH2_mMolar_DIFF': ['count', mean_diff_invlog(inc_low_n_vals=False)]}
        grouped_stats = self.test_foldchange_input_2.groupby(['FRAG_L', 'FRAG_R'], sort=False).agg(f)
        grouped_stats.columns = ['_'.join(col).strip() for col in grouped_stats.columns.values]
        
        grouped_stats = pd.DataFrame(grouped_stats).reset_index()
        grouped_stats.sort_index(axis=1, inplace=True)

        pd.util.testing.assert_frame_equal(grouped_stats, self.test_foldchange_goldenout_2)

    def test_mean_diff(self):

        f = {'Solubility_pH2_mMolar_DIFF': ['count', mean_diff(inc_low_n_vals=False)]}
        grouped_stats = self.test_foldchange_input_2.groupby(['FRAG_L', 'FRAG_R'], sort=False).agg(f)
        grouped_stats.columns = ['_'.join(col).strip() for col in grouped_stats.columns.values]

        grouped_stats = pd.DataFrame(grouped_stats).reset_index()
        grouped_stats.sort_index(axis=1, inplace=True)

        pd.util.testing.assert_frame_equal(grouped_stats, self.test_mean_diff_goldenout)

    def test_categorical_data(self):
        
        f = {'CAT_DATA_DIFF': [category_moved_up_one, category_moved_down_one, category_no_change,
                               category_changed_class, 'count']}
        grouped_stats = self.test_categorical_dataframe.groupby(['FRAG_L', 'FRAG_R'], sort=False).agg(f)
        grouped_stats.columns = ['_'.join(col).strip() for col in grouped_stats.columns.values]
        
        grouped_stats = pd.DataFrame(grouped_stats).reset_index()
        grouped_stats.sort_index(axis=1, inplace=True)

        # pp.pprint(grouped_stats)
        pd.util.testing.assert_frame_equal(grouped_stats, self.test_categorical_dataframe_golden)

    def test_categorical_data(self):

        f = {'CAT_DATA_DIFF': [category_moved_up_one_pct, category_moved_down_one_pct, category_no_change_pct,
                               category_changed_class_pct, 'count']}
        grouped_stats = self.test_categorical_dataframe.groupby(['FRAG_L', 'FRAG_R'], sort=False).agg(f)
        grouped_stats.columns = ['_'.join(col).strip() for col in grouped_stats.columns.values]

        grouped_stats = pd.DataFrame(grouped_stats).reset_index()
        grouped_stats.sort_index(axis=1, inplace=True)

        #pp.pprint(grouped_stats)
        pd.util.testing.assert_frame_equal(grouped_stats, self.test_categorical_dataframe_pct_golden)


if __name__ == '__main__':
    unittest.main()
