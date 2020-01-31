###################################################################
""" Summary: Class and Methods to manipulate MMP pairs

About: An object that reads a pairs file into memory for manipulation
by the associated methods.  Includes statistical aggregation functions
to summarise MMP data

"""
###################################################################
import logging
import sys
import os
import unittest
import tempfile
import pandas as pd
import numpy as np

import mmp.mmp_stats_functions as mmps
from mmp.mmp_data_objects import MMPDataObjectClass

def validate_agg_method(agg_method):
    """Validates that the keyword parameter passed to this module for the aggregation type
    is valid.  This can be imported and used by other scripts such as wrapper scripts. Because
    this might be called before the object is instantiated it is not an object level function."""

    agg_method = agg_method.upper()

    if agg_method == 'MEAN':
        return True

    # catch MEAN_DIFF before DIFFxx
    elif agg_method == 'MEAN_DIFF':
        return True

    elif agg_method == 'MEAN_DIFF_INVLOG':
        return True

    elif agg_method == 'CATEGORICAL':
        return True

    # order specific so should not get MEAN_DIFF
    elif 'DIFF' in agg_method:
        
        try:
            agg_method_int = int(agg_method.replace('DIFF', ''))
            if 100 > agg_method_int > 0:
                return True

            else:
                sys.exit('Invalid number for diffxx, try diffxx where 0 < xx < 100')

        except:
            sys.exit('Invalid parameter used for method diff, try diffxx where 0 < xx < 100')

    else:
        sys.exit('Invalid parameter used for method, try MEAN|MEAN_DIFF|MEAN_DIFF_INVLOG|CATEGORICAL|DIFFXX '
                 'where 0<XX<100')


def return_key_column_for_agg_method(column_name, agg_method):
    """Used by code like mmp_graph_predict to get the column name of the key aggregated data column
    after pandas based data aggregation.  For example, if I aggregate input column BLAH using MEAN_DIFF
    then I will expect the key return data column (delta for aggregated MMP's to be BLAH_MEAN_DIFF. This
    function is kept in same module as object mmp_pairs_object because it is better maintained alongside
    related function validate_agg_method and the actual object itself, keeping column renames consistent
    """

    agg_method = agg_method.upper()

    if agg_method == 'MEAN':
        return column_name + '_mean'

    elif agg_method == 'MEAN_DIFF':
        return column_name + '_mean_diff'

    elif agg_method == 'MEAN_DIFF_INVLOG':
        return column_name + '_mean_diff_involg'

    #
    # TODO: Need to work out how we use categorical data and therefore which column to use
    #
    elif agg_method == 'CATEGORICAL':
        return column_name + '_CATEGORICAL'

    # order specific so should not get MEAN_DIFF
    elif 'DIFF' in agg_method:

        try:
            agg_method_int = int(agg_method.replace('DIFF',''))
            if 100 > agg_method_int > 0:
                return column_name + '_diff' + str(agg_method_int)

            else:
                sys.exit('Invalid number for diffxx, try diffxx where 0 < xx < 100')

        except:
            sys.exit('Invalid parameter used for method diff, try diffxx where 0 < xx < 100')

    else:
        sys.exit('Invalid agg_method (%s) specified' % agg_method)


class MMPPairsObjectClass(MMPDataObjectClass):
    """Class implements objects nad methods for MMP pair manipulation"""

    def __init__(self, logger_object):
        """setup base object"""
        # core objects

        # core objects from base class
        MMPDataObjectClass.__init__(self, logger_object)

        # setup logger and check it's ok
        self.logger = logger_object
        if len(logging.Logger.manager.loggerDict) < 1:
            # exit with system status 1 and custom error
            sys.exit("Invalid or no logger object passed to MMPObjectClass.  Please create \
                    and pass a logger and set to use logging.disable if you don't want logging")
        
        # keep tabs on property data so we can aggregate MWT
        # and CLP using different stats than result data DIFFxx
        self.inc_props = False
        self.add_qmetric = False
        #
        self.grouped = {}

    def clear_out_data_pairobj(self):
        """Implements method from MMPDataObjectClass to clear out all data structures plus 
        clears out a few more related to this extended class"""

        self.clean_out_data()
        self.grouped.clear()

    def pairsdataobj_to_pd(self, cut_type, diff_col_name):
        """Iterate over the Pairs object and convert to a pandas dataframe for aggregation
        This is the only method in this class that uses methods and objects from the inherited 
        class MMPDataObjectsClass (mmp_data_objects)"""

        # fail if both single_pairs_dict and double_pairs_dict are empty
        if (len(self.single_pairs_dict) == 0) and (len(self.double_pairs_dict) == 0):
            self.logger.debug('No data found in single_pairs_dict and/or double_pairs_dict, expect no results')
            sys.exit("Error: no data found in single_pairs_dict and/or double_pairs_dict, nothing to find and write")

        # check cut_type, convert to int
        # TODO: code is repeated loads so should pull it out from everywhere and convert to function
        if cut_type.upper() == 'DOUBLE':
            # confusing but faster later
            cut_type_id = 3
        elif cut_type.upper() == 'BOTH':
            # confusing but faster later
            cut_type_id = 2
        elif cut_type.upper() == 'SINGLE':
            cut_type_id = 1
        else:
            self.logger.warn('cut_type specification is incorrect, using single cut: %s' % cut_type.upper())
            cut_type_id = 1

        # fail if diff_col_name has invalid value
        # print diff_col_name, self.headers_nosmi
        if diff_col_name not in self.headers_nosmi:
            sys.exit("Error: specified column name containing diff data does not exist")

        # find location of diff_col_name in the data table self.mol_data_dict
        #diff_col_name_idx = self.headers_nosmi.index(diff_col_name)

        headers = ['CUT', 'MOLID_L', 'MOLID_R', 'CONTEXT', 'FRAG_L', 'FRAG_R', 'ATTCHPT_CTX_L', 'ATTCHPT_FRAG_L',
                   'ATTCHPT_CTX_R', 'ATTCHPT_FRAG_R']

        for col in self.headers_nosmi:
            if col in list(self.headers_numeric_position.values()):
                headers.append(col)

        self.logger.info('Creating pandas dataframe from base mmp object via iterators')

        # Should be faster to pre-compute the dataframe size:
        # num_rows = 0
        # print "Final size: ", num_rows
        #
        # for el_ in self.single_pairs_dict.keys():
        #    size_ = len(self.single_pairs_dict[el_])
        #    num_rows = num_rows + size_
        # print "Final size: ", num_rows
        #
        # self.pairs_table = pd.DataFrame(index=np.arange(0, num_rows), columns=headers)
        # but hard due to pairs being dependent on factorial product of all pairs for each context,
        # less pairs generated via various conditions that get filtered out
        self.pairs_table = pd.DataFrame(columns=headers)

        row_position = 0

        # print pairs for single
        if cut_type_id <= 2:
            
            self.logger.debug("Adding single cuts to pandas dataframe")
            
            # first get pairs via iterator_single_pairs_dict...
            for molid_L, molid_R, ctx, frag_L, frag_R, fa_L, ca_L, fa_R, ca_R in self.iterator_single_pairs_dict():

                data_row = ['single', molid_L, molid_R, ctx, frag_L, frag_R, fa_L, ca_L, fa_R, ca_R]
    
                # ...now add the data differences
                # data_diff = self.mol_data_dict[molid_R][diff_col_name_idx] -
                #     self.mol_data_dict[molid_L][diff_col_name_idx]
                diff_array = self.get_data_diffs(molid_L, molid_R)
                data_row.extend(diff_array)
                
                # add row to table
                self.pairs_table.loc[row_position] = data_row
                row_position += 1

        # print pairs for double
        if cut_type_id >= 2:

            self.logger.debug("Adding double cuts to pandas dataframe")

            # first get pairs via iterator_double_pairs_dict...
            for molid_L, molid_R, ctx, frag_L, frag_R, fa_L, ca_L, fa_R, ca_R in self.iterator_double_pairs_dict():
                
                data_row = ['single', molid_L, molid_R, ctx, frag_L, frag_R, fa_L, ca_L, fa_R, ca_R]
                
                # ...now add the data differences
                # data_diff = self.mol_data_dict[molid_R][diff_col_name_idx] -
                #     self.mol_data_dict[molid_L][diff_col_name_idx]
                diff_array = self.get_data_diffs(molid_L, molid_R)
                # data_row.append(data_diff)
                data_row.extend(diff_array)

                # add row to table
                self.pairs_table.loc[row_position] = data_row
                row_position += 1

        self.logger.debug("Built dataframe of shape: %d, %d" % self.pairs_table.shape)

    def pd_read_csv(self, pairs_file, column_rename_dict = None):
        """Read in the CSV to a pandas dataframe:
        This method expects certain column names and if they do not exist in the input data 
        you need to rename them by providing a column_rename_dict on input parameter, see
        unit test for example by need columns:
        'FRAG_L', 'FRAG_R', 'ATTCHPT_FRAG_L', 'ATTCHPT_FRAG_R', 'DIFF'
        """

        # get and check CSV file
        self.csv_file = pairs_file
        if not os.path.isfile(self.csv_file):
            self.logger.warn('Cant instantiate object without a valid csv file input')
            sys.exit('Cant instantiate object without a valid csv file input')

        self.logger.info('Creating pandas dataframe from csv')
        if self.types_header:
            # Bug Fix:
            # In some cases we need to skip the second row in the CSV because we have a types header, this types header
            # is written only when -t flag invoked.  Read file back in with this types header
            # and the column types become object not float then aggregation fails.. use these lines before aggregation
            # to highlight the issue with an without -t flag: print self.pairs_table.dtypes && print self.grouped.dtypes
            # The self.types_header var had to be added to the class mmp_data_object
            self.pairs_table = pd.read_csv(self.csv_file, header=0, skiprows=[1], skipinitialspace=True,
                                           index_col=False)
        else:
            self.pairs_table = pd.read_csv(self.csv_file, header=0, skipinitialspace=True, index_col=False)

        # check to see if we need a column rename
        if column_rename_dict is not None:
            self.logger.info('Found a rename dict to use so will try to rename CSV columns')
            self.pairs_table.rename(columns = column_rename_dict, inplace=True)

        # check we have the columns we need later, these are hard coded later on so
        # must be dynamically renamed above or else everything will fail
        if 'MOLID_L' not in self.pairs_table.columns:
            self.logger.warn('Cannot find column MOLID_L in CSV file')
            sys.exit('Cannot find column MOLID_L in CSV file')

        elif 'MOLID_R' not in self.pairs_table.columns:
            self.logger.warn('Cannot find column MOLID_R in CSV file')
            sys.exit('Cannot find column MOLID_R in CSV file')
    
        elif 'FRAG_L' not in self.pairs_table.columns:
            self.logger.warn('Cannot find column FRAG_L in CSV file')
            sys.exit('Cannot find column FRAG_L in CSV file')

        elif 'FRAG_R' not in self.pairs_table.columns:
            self.logger.warn('Cannot find column FRAG_R in CSV file')
            sys.exit('Cannot find column FRAG_R in CSV file')

        elif 'ATTCHPT_FRAG_L' not in self.pairs_table.columns:
            self.logger.warn('Cannot find column ATTCHPT_FRAG_L in CSV file')
            sys.exit('Cannot find column ATTCHPT_FRAG_L in CSV file')

        elif 'ATTCHPT_FRAG_R' not in self.pairs_table.columns:
            self.logger.warn('Cannot find column ATTCHPT_FRAG_R in CSV file')
            sys.exit('Cannot find column ATTCHPT_FRAG_R in CSV file')

        else:
            self.logger.info('Column headers look great, proceed with aggregation')

        self.logger.info('Done reading from csv')

    def pd_aggregate_pairs_to_csv(self, filename, agg_type, agg_method='mean', prop_data=False,
                                  act_col=False, remove_id_dupes=True, inc_low_n_vals=True,
                                  add_qmetric=False):
        """Method to create grouped or aggregated pair data with associated stats and write to file.
        Can group data in different ways such as:
         frag = 'FRAG_L, FRAG_R' (to get the average change for a given fragment combination)
         attach = 'FRAG_L, FRAG_R, ATTCHPT_FRAG_L, ATTCHPT_FRAG_R' (to get the average change
             for for a given fragment and attachment point combination)
             
        Merged in a method (pairs_table.drop_duplicates) to remove ID/ID duplicates from a data frame.
        It's useful as a pre filter step before aggregation.

        The separate object level method pd_save_aggregate_pairs(self, out_file) got added back
        into this method due to memory issues (see inline comment below).  This plays off speed against 
        memory.  For small df's we can hold the aggregated pairs and .agg() data in memory then use .to_csv 
        on the whole df but for larger df's the script bombs with mem errors.  Thus for large df's we only
        generate the grouped object then iterate over it to calc the .agg() stats on an individual group at 
        a time, writing to file not mem. Small df = 100K rows and 15 columns, Large = 15,000K rows and 10 
        columns in testing so cutoff is an arbitary 100K * 15 = 1,500,000 (num_cells).  test_mode = True as
        input param will let you run the method designed for large df's over any small input df (testing mode)
        """

        # this is not fool proof.  If the incomming table does not have the columns ['MOL_L_CLP' 'MOL_R_CLP']
        # then code will pass here but fail later, so add second catch for missing columns
        if add_qmetric is True and prop_data is False:
            raise Exception('Must specify prop_data=True when requesting add_qmetric as Quality metric needs CLP data')
        # second catch - see above
        if add_qmetric is True:
            if 'MOL_L_CLP' not in self.pairs_table.columns or 'MOL_R_CLP' not in self.pairs_table.columns:
                raise Exception('Must have columns MOL_L_CLP and MOL_R_CLP present to add qmetric')

        def get_quality(row, clp_diff_idx, clp_diff_sd_idx, clp_left_sd_idx, num_pairs_idx):
            """ Function that takes list input (row) and returns an estimate of quality based on
             Good, Number of pairs >= 15 & stdev of L.clogp >= 0.5 & when |d.clogp|>0, stdev of d.clogp>0
             Medium, Number of pairs >= 6 and stdev of L.clogp >= 0.3 and when |d.clogp|>0, stdev of d.clogp>0
             Poor, The rest of the matched-pairs that do not satisfied the above criteria
             Columns to work with are
                'MOL_CLP_DIFF, 'MOL_CLP_DIFF_std', 'MOL_L_CLP_std', act_col + '_count'
            """
            # first write of function with simplest slower logic
            # if row[num_pairs_idx] >= 15 and row[clp_left_sd_idx] >= 0.5 and \
            #   abs(row[clp_diff_idx]) > 0 and row[clp_diff_sd_idx] > 0:
            #    return 2
            # elif row[num_pairs_idx] >= 6 and row[clp_left_sd_idx] >= 0.3 and \
            #   abs(row[clp_diff_sd_idx]) > 0 and row[clp_diff_sd_idx] > 0:
            #    return 1
            # else:
            #    return 0
            if abs(row[clp_diff_idx]) > 0 and row[clp_diff_sd_idx] > 0:

                num_pairs = row[num_pairs_idx]
                clp_left_sd = row[clp_left_sd_idx]

                if num_pairs >= 15 and clp_left_sd >= 0.5:
                    return 2

                elif num_pairs >= 6 and clp_left_sd >= 0.3:
                    return 1

                else:
                    return 0
            else:
                return 0

        def write_grouped_table_to_csv(filename):
            # This method is a work around a bug seen with multiline headers in pd 0.16.1 with py 2.7 I cannot write a
            # multi line header: https://github.com/pydata/pandas/issues/5539 because I get an additional empty header
            # line in the csv.  Instead we write the original header to csv, replace it then append the whole df to csv

            # 201701-26 Adding Quality metric
            # TODO: Putting this aggregate function here in a submethod is messy but simplest approach for now
            #       execution suffers badly and it needs a refactor
            if add_qmetric:
                self.logger.info('Calculating Quality metric')
                header_cols = list(self.grouped.columns)
                clp_diff_idx = header_cols.index('MOL_CLP_DIFF_mean')
                clp_diff_sd_idx = header_cols.index('MOL_CLP_DIFF_std')
                clp_left_sd_idx = header_cols.index('MOL_L_CLP_std')
                num_pairs_idx = header_cols.index(act_col + '_count')
                self.grouped['QUALITY'] = self.grouped.apply(lambda row: get_quality(row, clp_diff_idx,
                                                                                     clp_diff_sd_idx, clp_left_sd_idx,
                                                                                     num_pairs_idx), axis=1)

            if self.types_header:
                self.logger.info('Writing Types Header')
                # write the header to file
                index_count = len(self.grouped.index.names)
                for idx in range(0, index_count):
                    self.grouped.reset_index(level=0, inplace=True)
                pd.DataFrame(data=[self.grouped.columns]).to_csv(filename,
                                                                 float_format='%.3f',
                                                                 header=False,
                                                                 index=False,
                                                                 na_rep="NaN")

                # now replace header with column type labels
                types_header_for_insert = list(self.grouped.columns.values)
                for idx, val in enumerate(self.grouped.columns.values):
                    # print idx, val, self.grouped[val].dtype
                    if self.grouped[val].dtype == 'float64':
                        types_header_for_insert[idx] = 'REAL'

                    elif self.grouped[val].dtype == 'int64':
                        types_header_for_insert[idx] = 'INTEGER'

                    else:
                        types_header_for_insert[idx] = 'STRING'

                self.grouped.columns = types_header_for_insert

                # and append the whole df with new header to original file
                self.logger.info('Writing to CSV')
                self.grouped.to_csv(filename, mode="a", float_format='%.3f', index=False, na_rep="NaN")

            else:
                # simple write
                self.logger.info('Writing to CSV')
                self.grouped.to_csv(filename, float_format='%.3f', na_rep="NaN")

        # As the DIFF and other functions can be expensive, we will only allow
        agg_method = agg_method.upper()
        agg_type = agg_type.upper()
        out_file = filename
        tmp_out_file = out_file + "_tmp"
        
        validate_agg_method(agg_method)

        # will be 0 if 'DIFFXX' or >0 if 'MEAN_DIFF' and -1 when string does not contain DIFF
        diff_type = str.find(agg_method, 'DIFF')

        if act_col:
            act_col += '_DIFF'

        # As the DIFF and other functions can be expensive, we will only allow it to be used with act_col
        # test times gave 22 sec for mean versus 202 sec for diffxx over 1.4M keys, single column
        if act_col is False and 'DIFF' in agg_method:
            sys.exit('Please specify the activity difference column name when using the DIFFXX function as it is '
                     'expensive to run')

        # drop columns we don't need:
        self.logger.info('Drop columns we dont need')

        if prop_data:
            self.pairs_table.rename(columns={'W_AMW_DIFF': 'MOL_MWT_DIFF'}, inplace=True)

        if act_col:
            if prop_data:
                column_list = ['MOLID_L', 'MOLID_R', 'FRAG_L', 'FRAG_R', 'ATTCHPT_CTX_L',
                               'ATTCHPT_CTX_R', 'MOL_CLP_DIFF', 'MOL_MWT_DIFF', act_col]
                if add_qmetric:
                    column_list.extend(['MOL_L_CLP', 'MOL_R_CLP'])
                self.pairs_table = self.pairs_table[column_list]

            else:
                column_list = ['MOLID_L', 'MOLID_R', 'FRAG_L', 'FRAG_R', 'ATTCHPT_CTX_L', 'ATTCHPT_CTX_R', act_col]
                self.pairs_table = self.pairs_table[column_list]

        # could sort table before we remove duplicates to get consistent smiles ouptut
        # tested on 75K input SMI but cost was high
        # self.logger.info('sort data table')
        # self.pairs_table.sort(['MOLID_L', 'MOLID_R', 'FRAG_L', 'FRAG_R', 'ATTCHPT_CTX_L', 'ATTCHPT_CTX_R'],
        #    inplace=True)
        
        # drop duplicates
        self.logger.info('Drop duplicate rows')
        
        if agg_type == 'FRAG':
            if remove_id_dupes:
                self.pairs_table.drop_duplicates(['MOLID_L', 'MOLID_R', 'FRAG_L', 'FRAG_R'], inplace=True)
            self.pairs_table.drop(['MOLID_L', 'MOLID_R', 'ATTCHPT_CTX_L', 'ATTCHPT_CTX_R'], axis=1, inplace=True)

        elif agg_type == 'ATTACH':
            # 
            # This creates a situation where a different N (count) can occur for non attachment point aggregated data
            # versus attachment point aggregated.
            if remove_id_dupes:
                self.pairs_table.drop_duplicates(['MOLID_L', 'MOLID_R', 'FRAG_L', 'FRAG_R',
                                                 'ATTCHPT_CTX_L', 'ATTCHPT_CTX_R'],
                                                 inplace=True)
            # as opposed to this, which will select an arbitary one of the many attachment points for the given pairs
            # self.pairs_table.drop_duplicates(['MOLID_L', 'MOLID_R', 'FRAG_L', 'FRAG_R'], inplace = True)
            self.pairs_table.drop(['MOLID_L', 'MOLID_R'], axis=1, inplace=True)

        else:
            sys.exit('Invalid value for parameter agg_type: %s (should be FRAG or ATTACH' % agg_type)

        # do grouping and aggregation
        self.logger.info('Group')

        # python3 needs numeric not object for aggregate
        self.pairs_table = self.pairs_table.apply(pd.to_numeric, errors='ignore')

        if agg_type == 'ATTACH':
            self.grouped = self.pairs_table.groupby(['FRAG_L', 'FRAG_R', 'ATTCHPT_CTX_L', 'ATTCHPT_CTX_R'], sort=False)

        # assume we have caught all miss-assignments of agg_type so must have 'FRAG'
        else:
            self.grouped = self.pairs_table.groupby(['FRAG_L', 'FRAG_R'], sort=False)

        ##########################
        #
        # See method level comments regarding df size. pd.df's of 15,000K rows failed (~80K input SMI)
        # Therefore two method implemented to balance speed versus memory
        # TODO: must be a better way to deal with below multi switch possibilities...
        #
        #############################
        # (height, width) = self.pairs_table.shape
        # num_cells = height * width

        self.logger.info('Aggregating in memory')
        if agg_method == 'MEAN':
            
            self.logger.info('Aggregating using MEAN')
            self.grouped = self.grouped.aggregate(['mean', 'std', 'count'])
            self.grouped.columns = ['_'.join(col).strip() for col in self.grouped.columns.values]

            # try:
            write_grouped_table_to_csv(filename)
            # except:
            #    self.logger.warn('You have not done any aggregation, use aggregate_pairs() first')
            #    sys.exit('You have not done any aggregation, use aggregate_pairs() first')

        elif agg_method == 'CATEGORICAL':
            
            self.logger.info('Aggregating CATEGORICAL data')

            if prop_data:
                f = {act_col: ['count', mmps.category_moved_up_one_pct, mmps.category_moved_down_one_pct,
                               mmps.category_no_change_pct, mmps.category_changed_class_pct],
                     'MOL_CLP_DIFF': ['mean', 'std'],
                     'MOL_MWT_DIFF': ['first']}
                if add_qmetric:
                    f['MOL_L_CLP'] = ['std']

            else:
                f = {act_col: ['count', mmps.category_moved_up_one_pct, mmps.category_moved_down_one_pct,
                               mmps.category_no_change_pct, mmps.category_changed_class_pct]}

            self.grouped = self.grouped.aggregate(f)
            self.grouped.columns = ['_'.join(col).strip() for col in self.grouped.columns.values]

            write_grouped_table_to_csv(filename)

        # case where we have MEAN_DIFF or DIFFxx
        elif diff_type >= 0:

            if diff_type == 0:
                try:
                    agg_method_int = int(agg_method.replace('DIFF', ''))

                except:
                    sys.exit('Invalid parameter used for method diff, try diffxx where 0 < xx < 100')
            
            if prop_data:
                
                self.logger.info('Aggregating using DIFF or MEAN_DIFF, with property data')
                # specify what to use on what column
                if diff_type == 0:
                    f = {act_col: [mmps.diffn_list_rtn(agg_method_int, inc_low_n_vals=inc_low_n_vals), 'count',
                                   mmps.n_pos_diff, mmps.n_neg_diff], 'MOL_CLP_DIFF': ['mean', 'std'],
                         'MOL_MWT_DIFF': ['first']}
                    if add_qmetric:
                        f['MOL_L_CLP'] = ['std']
                else:
                    if 'INVLOG' in agg_method:
                        f = {act_col: [mmps.mean_diff_invlog(inc_low_n_vals=inc_low_n_vals), 'count', mmps.n_pos_diff,
                                       mmps.n_neg_diff], 'MOL_CLP_DIFF': ['mean', 'std'],
                             'MOL_MWT_DIFF': ['first']}
                        if add_qmetric:
                            f['MOL_L_CLP'] = ['std']
                    else:
                        f = {act_col: [mmps.mean_diff(inc_low_n_vals=inc_low_n_vals), 'count', mmps.n_pos_diff,
                                       mmps.n_neg_diff], 'MOL_CLP_DIFF': ['mean', 'std'],
                             'MOL_MWT_DIFF': ['first']}
                        if add_qmetric:
                            f['MOL_L_CLP'] = ['std']

                self.grouped = self.grouped.aggregate(f)
                self.grouped.columns = ['_'.join(col).strip() for col in self.grouped.columns.values]
                # self.grouped[[act_col+'_'+agg_method+'_low', act_col+'_'+agg_method,
                # act_col+'_'+agg_method+'_upp']] = self.grouped[act_col+'_'+agg_method.lower()+'_all'].apply(pd.Series)
                # self.grouped.drop(act_col+'_'+agg_method.lower()+'_all', axis=1, inplace=True)

                # self.grouped.dropna(how='any', inplace=True)
                write_grouped_table_to_csv(tmp_out_file)

            else:

                self.logger.info('Aggregating using DIFF, with no property data')

                if diff_type == 0:
                    f = {act_col: ['count', mmps.n_pos_diff, mmps.n_neg_diff,
                                   mmps.diffn_list_rtn(agg_method_int, inc_low_n_vals=inc_low_n_vals)]}
                else:
                    if 'INVLOG' in agg_method:
                        f = {act_col: ['count', mmps.n_pos_diff, mmps.n_neg_diff,
                                       mmps.mean_diff_invlog(inc_low_n_vals=inc_low_n_vals)]}
                    else:
                        f = {act_col: ['count', mmps.n_pos_diff, mmps.n_neg_diff,
                                       mmps.mean_diff(inc_low_n_vals=inc_low_n_vals)]}

                self.grouped = self.grouped.aggregate(f)

                self.grouped.columns = ['_'.join(col).strip() for col in self.grouped.columns.values]

                write_grouped_table_to_csv(tmp_out_file)

            # could not get efficient way to break columns out in memory or via group iterator so doing it the ugly way
            # need to replace ,ACT_A_DIFF_diff60_all, with ,ACT_A_DIFF_DIFF60_low,ACT_A_DIFF_DIFF60,ACT_A_DIFF_DIFF60_upp
            # and sort out ,"(39.999999751010833, 39.99999997345952, 39.999999997170974)",
            # TODO: Get rid of this file based cleanup and implement method for DIFFXX calc to return separate columns
            #
            self.logger.info('Messy cleanup of function output')

            if diff_type == 0:
                header_string_orig = act_col + '_' + agg_method.lower() + '_all'
            else:
                header_string_orig = act_col + '_' + agg_method.lower()

            header_string_new = act_col + '_' + agg_method.lower() + '_low' + ', ' + act_col + '_' + \
                                agg_method.lower() + ', '  + act_col + '_' + agg_method.lower() + '_upp'

            with open(tmp_out_file, 'rt') as f_read:
                with open(out_file, 'wt') as f_write:
                    line_num = 0
                    for x in f_read:
                        line_num += 1
                        y = x
                        if line_num == 1:
                            y = y.replace(header_string_orig, header_string_new)
                            y = y.replace(" ", '')
                        elif line_num == 2 and self.types_header:
                            # strip eol
                            y = y.rstrip("\n")
                            # to list then find all 'STRING', check the underlying data, if (nan, nan, nan), change to
                            #  REAL
                            x = y.split(",")
                            for idx, item in enumerate(x):
                                if item == 'STRING':
                                    # print "->", idx, item, self.grouped.iloc[2, idx]
                                    # does this column position correspond to the (nan, nan, nan) position in dataframe
                                    if str(self.grouped.iloc[2, idx]).find(",") > 0:
                                        # now replace with REAL,REAL,REAL
                                        # print "--->", idx, item, self.grouped.iloc[2, idx]
                                        x[idx] = 'REAL,REAL,REAL'
                            # convert back to string
                            y = ','.join(x)
                            y += "\n"
                        else:
                            y = x
                            y = y.replace(",\"[", ',')
                            y = y.replace("]\",", ',')
                            y = y.replace("NaN]\"", 'NaN')
                            y = y.replace(",\"(", ',')
                            y = y.replace(")\",", ',')
                            y = y.replace(")\"", '')
                            y = y.replace(", ", ',')
                            # y = y.replace(",,",',NaN,')
                            # y = y.replace(", ,",',NaN,')
                        # Java needs NaN not nan - not sure what we'll do if we remove this loop for native to_csv
                        y = y.replace(",nan", ',NaN')
                        f_write.write(y)

            os.remove(tmp_out_file)
        
        else:
            sys.exit("Invalid input for agg_method")

        self.logger.info('Completed aggregation')     


#
# unittest everything
#
class _Test_MMPPairsObjectClass(unittest.TestCase):
    """Test class to test the object and methods"""

    def setUp(self):
        #
        self.maxDiff = None

        # setup test data location use tempfile.NamedTemporaryFile(delete=False) to persist data on disk
        self.temp_file_input_pairs = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='w+t')
        self.temp_file_input_pairdupes = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='w+t')
        self.temp_file_input_diffdata = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='w+t')
        self.temp_file_input_diffdata_02 = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='w+t')
        self.temp_file_output = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='w+t')
        self.temp_file_output_nodupes = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='w+t')
        self.temp_file_output_diffdata = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='w+t')
        self.temp_file_input_csv = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='w+t')

        self.mmplogger = logging.getLogger('mmpobjectclass_testlogger')
        logging.disable(logging.CRITICAL)

        self.test_mmp_pairs_object = MMPPairsObjectClass(self.mmplogger)
        # TODO: move files write blocks into each unittest

        self.test_dataset_goldeninput_header = 'MOLID_L, MOLID_R, CONTEXT, FRAG_L, FRAG_R, ATTCHPT_CTX_L, '
        self.test_dataset_goldeninput_header += 'ATTCHPT_FRAG_L, ATTCHPT_CTX_R, ATTCHPT_FRAG_R, pIC50_DIFF'

        # I need to remap my input data columns to something consistent
        self.test_columns_dict = {'FRAG_L' : 'FRAG_L',
                                  'FRAG_R' : 'FRAG_R',
                                  'ATTCHPT_FRAG_L' : 'ATTCHPT_FRAG_L',
                                  'ATTCHPT_FRAG_R' : 'ATTCHPT_FRAG_R',
                                  'pIC50_DIFF' : 'DIFF'
                                 }

        self.test_dataframe = pd.DataFrame({
            'A' : ['aa', 'ab', 'ac', 'ad', 'aa', 'ab', 'ac', 'ad', 'aa', 'aa', 'ab'],
            'B' : [-0.542, -3.043, 0.264, 0.094, -0.262, 0.344, 0.094, -0.262, -0.555, -0.54, -0.27],
            'C' : [-3.043, 0.264, 0.094, -0.262, 0.344, 0.769, 0.094, -0.262, -3.001, -3.10, 0.35],
            'D' : [0.264, 0.094, -0.262, 0.344, 0.769, 0.811, 0.094, -0.262, -0.260, -3.001, 0.100],
            'E' : [-0.262, 0.344, 0.769, 0.811, -1.350, -1.475, 0.094, -0.262, -0.254, -0.254, 0.901],
            'F' : [-0.262, 0.344, 0.769, 0.811, -1.350, -1.475, 0.094, -0.262, -0.256, -0.206, 0.344]
        })
        
        self.test_dataframe_result = pd.DataFrame({
            'A': ['aa', 'ab', 'ac', 'ad'],
            'diff60': [-59.06244, 21.25944, np.nan, np.nan],
            'len': [4.0, 3.0, 2.0, 2.0],
            'mean': [-2.200, 0.461, 0.094, -0.262]
        })

        self.test_dataset_goldeninput_csv_headers = \
                ['SMILES','PIC50','CHEMBL_ID','CHEMBL_TARGET_ID','CHEMBL_TARGET']

        self.test_dataset_goldeninput_csv_data = {
            # Below data is taken from DOI: 10.12688/f1000research.3-36.v2
            # -> https://f1000research.com/articles/3-36
            # Bajorath et. al. 2014 "Matched molecular pair-based data sets for computer-aided medicinal chemistry."
            # All ID's are CHEMBL Compound or Target IDs
            'O(CC)c1ccc(N2C(=Nc3c(cccc3)C2=O)[C@H](N(C(=O)Cc2ccc(cc2)-c2ccccc2)CCOCC)CC)cc1,8.045757491,230126,4441,C-X-C chemokine receptor type 3': None,
            'O(CC)c1ccc(N2C(=Nc3c(cccc3)C2=O)[C@H](N(C(=O)Cc2ccc(cc2)-c2ccccc2)CCOCC)C)cc1,7.124938737,266663,4441,C-X-C chemokine receptor type 3': None,
            'O(CC)c1ccc(N2C(=Nc3c(cccc3)C2=O)CN(C(=O)Cc2ccc(cc2)-c2ccccc2)CCOCC)cc1,5.638272164,231487,4441,C-X-C chemokine receptor type 3': None,
            'O(CC)c1ccc(N2C(=Nc3c(cccc3)C2=O)[C@H](N(C(=O)Cc2ccc(cc2)-c2ccccc2)CCOCC)c2ccccc2)cc1,5.397940009,230127,4441,C-X-C chemokine receptor type 3': None,
            '[nH]1nc(c(c1)-c1cc(ncc1)-c1ccc(cc1)C)-c1ncccc1,7.552841969,205156,4439,TGF-beta receptor type I': None,
            '[nH]1nc(c(c1)-c1cc(ncc1)-c1cc(ccc1)C)-c1ncccc1,6.718966633,205260,4439,TGF-beta receptor type I': None,
            '[nH]1nc(c(c1)-c1cc(ncc1)-c1ccccc1C)-c1ncccc1,5.400007822,383433,4439,TGF-beta receptor type I': None,
            'Clc1cc(cnc1N1C[C@@H]([NH+](CC1)C1CC[NH+](CC1)Cc1ccc(Cl)cc1)CC)C(=O)NC,8.522878745,1681874,4441,C-X-C chemokine receptor type 3': None,
            'Clc1cc(cnc1N1C[C@@H]([NH+](CC1)C1CC[NH+](CC1)Cc1ccc(Cl)cc1)C)C(=O)NC,7.494850022,1681841,4441,C-X-C chemokine receptor type 3': None,
            'Clc1cc(cnc1N1CC[NH+](CC1)C1CC[NH+](CC1)Cc1ccc(Cl)cc1)C(=O)NC,6.251811973,565761,4441,C-X-C chemokine receptor type 3': None,
            'Clc1cc(cnc1N1C[C@@H]([NH+](CC1)C1CC[NH+](CC1)Cc1ccc(Cl)cc1)c1ccccc1)C(=O)NC,6.080921908,1681876,4441,C-X-C chemokine receptor type 3': None,
            'FC(F)(F)c1ccc(cc1)-c1nccc(c1)-c1c[nH]nc1-c1ncccc1,7.619788758,383523,4439,TGF-beta receptor type I': None,
            'FC(F)(F)c1cc(ccc1)-c1nccc(c1)-c1c[nH]nc1-c1ncccc1,6.619788758,382466,4439,TGF-beta receptor type I': None,
            'FC(F)(F)c1ccccc1-c1nccc(c1)-c1c[nH]nc1-c1ncccc1,5.029978734,426852,4439,TGF-beta receptor type I': None
        }

        # example pairs file (CHEMBL data):
        # below output is derived from command line given, with testdata.csv file contents from above TGF-beta lines:
        # getMMPStatsfromCSV.sh -i testdata.csv -o testdata.pairs -s SMILES -n CHEMBL_ID -a PIC50 -c SINGLE
        self.test_dataset_goldeninput_pairs = {
'single,205156,205260,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1cc[1cH]cc1,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-8.339e-01': None,
'single,205156,383433,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1cc[1cH]cc1,Cc1[1cH]cccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-2.153e+00': None,
'single,205156,383523,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1cc[1cH]cc1,FC(c1cc[1cH]cc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],6.695e-02': None,
'single,205156,382466,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1cc[1cH]cc1,FC(c1c[1cH]ccc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-9.331e-01': None,
'single,205156,426852,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1cc[1cH]cc1,FC(c1[1cH]cccc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-2.523e+00': None,
'single,205260,205156,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1c[1cH]ccc1,Cc1cc[1cH]cc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],8.339e-01': None,
'single,205260,383433,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1c[1cH]ccc1,Cc1[1cH]cccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-1.319e+00': None,
'single,205260,383523,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1c[1cH]ccc1,FC(c1cc[1cH]cc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],9.008e-01': None,
'single,205260,382466,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1c[1cH]ccc1,FC(c1c[1cH]ccc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-9.918e-02': None,
'single,205260,426852,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1c[1cH]ccc1,FC(c1[1cH]cccc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-1.689e+00': None,
'single,383433,205156,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1[1cH]cccc1,Cc1cc[1cH]cc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],2.153e+00': None,
'single,383433,205260,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1[1cH]cccc1,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],1.319e+00': None,
'single,383433,383523,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1[1cH]cccc1,FC(c1cc[1cH]cc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],2.220e+00': None,
'single,383433,382466,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1[1cH]cccc1,FC(c1c[1cH]ccc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],1.220e+00': None,
'single,383523,205156,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1cc[1cH]cc1)(F)F,Cc1cc[1cH]cc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-6.695e-02': None,
'single,383523,205260,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1cc[1cH]cc1)(F)F,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-9.008e-01': None,
'single,383523,383433,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1cc[1cH]cc1)(F)F,Cc1[1cH]cccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-2.220e+00': None,
'single,383523,382466,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1cc[1cH]cc1)(F)F,FC(c1c[1cH]ccc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-1.000e+00': None,
'single,383523,426852,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1cc[1cH]cc1)(F)F,FC(c1[1cH]cccc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-2.590e+00': None,
'single,382466,205156,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1c[1cH]ccc1)(F)F,Cc1cc[1cH]cc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],9.331e-01': None,
'single,382466,205260,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1c[1cH]ccc1)(F)F,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],9.918e-02': None,
'single,382466,383523,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1c[1cH]ccc1)(F)F,FC(c1cc[1cH]cc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],1.000e+00': None,
'single,382466,426852,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1c[1cH]ccc1)(F)F,FC(c1[1cH]cccc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-1.590e+00': None,
'single,426852,205156,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1[1cH]cccc1)(F)F,Cc1cc[1cH]cc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],2.523e+00': None,
'single,426852,205260,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1[1cH]cccc1)(F)F,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],1.689e+00': None,
'single,426852,383523,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1[1cH]cccc1)(F)F,FC(c1cc[1cH]cc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],2.590e+00': None,
'single,426852,382466,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1[1cH]cccc1)(F)F,FC(c1c[1cH]ccc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],1.590e+00': None,
'single,205156,383523,[nH]1[n]c(c(c1)c1cc([n]cc1)c1cc[1cH]cc1)c1[n]cccc1,[1CH4],F[1CH](F)F,[1:CAR],[1:C3],[1:CAR],[1:C3],6.695e-02': None,
'single,383523,205156,[nH]1[n]c(c(c1)c1cc([n]cc1)c1cc[1cH]cc1)c1[n]cccc1,F[1CH](F)F,[1CH4],[1:CAR],[1:C3],[1:CAR],[1:C3],-6.695e-02': None,
'single,205260,382466,[nH]1[n]c(c(c1)c1cc([n]cc1)c1c[1cH]ccc1)c1[n]cccc1,[1CH4],F[1CH](F)F,[1:CAR],[1:C3],[1:CAR],[1:C3],-9.918e-02': None,
'single,382466,205260,[nH]1[n]c(c(c1)c1cc([n]cc1)c1c[1cH]ccc1)c1[n]cccc1,F[1CH](F)F,[1CH4],[1:CAR],[1:C3],[1:CAR],[1:C3],9.918e-02': None,
'single,382466,383433,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1c[1cH]ccc1)(F)F,Cc1[1cH]cccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-1.220e+00': None,
'single,383433,426852,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1[1cH]cccc1,FC(c1[1cH]cccc1)(F)F,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-3.700e-01': None,
'single,426852,383433,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1[1cH]cccc1)(F)F,Cc1[1cH]cccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],3.700e-01': None,
'single,383433,426852,[nH]1[n]c(c(c1)c1cc(c2[1cH]cccc2)[n]cc1)c1[n]cccc1,[1CH4],F[1CH](F)F,[1:CAR],[1:C3],[1:CAR],[1:C3],-3.700e-01': None,
'single,426852,383433,[nH]1[n]c(c(c1)c1cc(c2[1cH]cccc2)[n]cc1)c1[n]cccc1,F[1CH](F)F,[1CH4],[1:CAR],[1:C3],[1:CAR],[1:C3],3.700e-01': None
        }

        # These two CHEMBL mols will produce very repetitive pairs. On aggregation we should only count
        # the id-id pair once if filtering is on, see test
        # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL268389
        # CC\C(=C(/c1ccc(I)cc1)\c2ccc(OCCCCCN3CCCC3)cc2)\c4ccccc4,6.1,268389,2095150,Phosphodiesterase 1
        # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL267035/
        # CC\C(=C(/c1ccc(I)cc1)\c2ccc(OCCCCN3CCCC3)cc2)\c4ccccc4,6.0,267035,2095150,Phosphodiesterase 1
        #        ['SMILES','PIC50','CHEMBL_ID','CHEMBL_TARGET_ID','CHEMBL_TARGET']
        # a lot of data is removed as unnecessary to trigger test condition

        self.test_dataset_goldeninput_pairdupes_header = 'CUT,MOLID_L,MOLID_R,CONTEXT,FRAG_L,FRAG_R,ATTCHPT_CTX_L,ATTCHPT_FRAG_L,ATTCHPT_CTX_R,ATTCHPT_FRAG_R,PIC50_DIFF'

        self.test_dataset_goldeninput_pairdupes = {
            'double,268389,267035,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc([1OH])cc2)cc1.[2CH3]CCN1CCCC1,[2CH3][1CH3],[12CH4],[2:C3|1:O3],[1:C3|2:C3],[2:C3|1:O3],[1:C3|2:C3],-1.000e-01': None,
            'double,268389,267035,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc([2OH])cc2)cc1.[1CH3]CCN1CCCC1,[1CH3][2CH3],[12CH4],[1:C3|2:O3],[2:C3|1:C3],[1:C3|2:O3],[2:C3|1:C3],-1.000e-01': None,
            'double,267035,268389,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc([1OH])cc2)cc1.[2CH3]CCN1CCCC1,[12CH4],[2CH3][1CH3],[2:C3|1:O3],[1:C3|2:C3],[2:C3|1:O3],[1:C3|2:C3],1.000e-01': None,
            'double,267035,268389,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc([2OH])cc2)cc1.[1CH3]CCN1CCCC1,[12CH4],[1CH3][2CH3],[1:C3|2:O3],[2:C3|1:C3],[1:C3|2:O3],[2:C3|1:C3],1.000e-01': None,
            'double,268389,267035,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(O[1CH3])cc2)cc1.[2CH3]CN1CCCC1,[2CH3][1CH3],[12CH4],[2:C3|1:C3],[1:C3|2:C3],[2:C3|1:C3],[1:C3|2:C3],-1.000e-01': None,
            'double,268389,267035,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(O[2CH3])cc2)cc1.[1CH3]CN1CCCC1,[1CH3][2CH3],[12CH4],[1:C3|2:C3],[2:C3|1:C3],[1:C3|2:C3],[2:C3|1:C3],-1.000e-01': None,
            'double,267035,268389,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(O[1CH3])cc2)cc1.[2CH3]CN1CCCC1,[12CH4],[2CH3][1CH3],[2:C3|1:C3],[1:C3|2:C3],[2:C3|1:C3],[1:C3|2:C3],1.000e-01': None,
            'double,267035,268389,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(O[2CH3])cc2)cc1.[1CH3]CN1CCCC1,[12CH4],[1CH3][2CH3],[1:C3|2:C3],[2:C3|1:C3],[1:C3|2:C3],[2:C3|1:C3],1.000e-01': None,
            'double,268389,267035,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(O[1CH3])cc2)cc1.[2CH3]N1CCCC1,[2CH3]C[1CH3],[2CH3][1CH3],[2:C3|1:C3],[1:C3|2:C3],[2:C3|1:C3],[1:C3|2:C3],-1.000e-01': None,
            'double,268389,267035,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(O[2CH3])cc2)cc1.[1CH3]N1CCCC1,[1CH3]C[2CH3],[1CH3][2CH3],[1:C3|2:C3],[2:C3|1:C3],[1:C3|2:C3],[2:C3|1:C3],-1.000e-01': None,
            'double,267035,268389,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(O[1CH3])cc2)cc1.[2CH3]N1CCCC1,[2CH3][1CH3],[2CH3]C[1CH3],[2:C3|1:C3],[1:C3|2:C3],[2:C3|1:C3],[1:C3|2:C3],1.000e-01': None,
            'double,267035,268389,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(O[2CH3])cc2)cc1.[1CH3]N1CCCC1,[1CH3][2CH3],[1CH3]C[2CH3],[1:C3|2:C3],[2:C3|1:C3],[1:C3|2:C3],[2:C3|1:C3],1.000e-01': None,
            'double,268389,267035,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(OC[1CH3])cc2)cc1.[2CH3]N1CCCC1,[2CH3][1CH3],[12CH4],[2:C3|1:C3],[1:C3|2:C3],[2:C3|1:C3],[1:C3|2:C3],-1.000e-01': None,
            'double,268389,267035,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(OC[2CH3])cc2)cc1.[1CH3]N1CCCC1,[1CH3][2CH3],[12CH4],[1:C3|2:C3],[2:C3|1:C3],[1:C3|2:C3],[2:C3|1:C3],-1.000e-01': None,
            'double,267035,268389,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(OC[1CH3])cc2)cc1.[2CH3]N1CCCC1,[12CH4],[2CH3][1CH3],[2:C3|1:C3],[1:C3|2:C3],[2:C3|1:C3],[1:C3|2:C3],1.000e-01': None,
            'double,267035,268389,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(OC[2CH3])cc2)cc1.[1CH3]N1CCCC1,[12CH4],[1CH3][2CH3],[1:C3|2:C3],[2:C3|1:C3],[1:C3|2:C3],[2:C3|1:C3],1.000e-01': None,
            'double,268389,267035,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(OCC[1CH3])cc2)cc1.[2NH]1CCCC1,[2CH3][1CH3],[12CH4],[2:N3|1:C3],[1:C3|2:C3],[2:N3|1:C3],[1:C3|2:C3],-1.000e-01': None,
            'double,268389,267035,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(OCC[2CH3])cc2)cc1.[1NH]1CCCC1,[1CH3][2CH3],[12CH4],[1:N3|2:C3],[2:C3|1:C3],[1:N3|2:C3],[2:C3|1:C3],-1.000e-01': None,
            'double,267035,268389,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(OCC[1CH3])cc2)cc1.[2NH]1CCCC1,[12CH4],[2CH3][1CH3],[2:N3|1:C3],[1:C3|2:C3],[2:N3|1:C3],[1:C3|2:C3],1.000e-01': None,
            'double,267035,268389,Ic1ccc(C(=C(CC)c2ccccc2)c2ccc(OCC[2CH3])cc2)cc1.[1NH]1CCCC1,[12CH4],[1CH3][2CH3],[1:N3|2:C3],[2:C3|1:C3],[1:N3|2:C3],[2:C3|1:C3],1.000e-01': None
        }

        #
        # fake test data, numeric id's for compounds/fragments have not relation to real data
        self.test_dataset_goldeninput_diffheader = 'MOLID_L, MOLID_R, FRAG_L, FRAG_R, ATTCHPT_CTX_L, ATTCHPT_FRAG_L, ATTCHPT_CTX_R, ATTCHPT_FRAG_R, DIFF_DIFF'

        self.test_dataset_goldeninput_diffdata = {
'19733, 19733, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -0.542,': None,
'88154, 94170, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -3.043,': None,
'22633, 20295, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], 0.264,': None,
'07788, 08310, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], 0.094,': None,
'00035, 03079, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -0.262,': None,
'82657, 87490, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], 0.344,': None,
'82472, 87490, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], 0.769,': None,
'82657, 87490, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], 0.811,': None,
'04116, 04139, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -1.35,': None,
'05358, 05056, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -1.475,': None,
'32144, 32058, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -0.027,': None,
'93314, 05098, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -3.2,': None,
'91240, 05098, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -3.345,': None,
'93399, 05098, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.95,': None,
'16563, 16207, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -0.107,': None,
'16278, 16060, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.298,': None,
'88346, 94030, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -3.273,': None,
'18283, 20381, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.301,': None,
'21549, 21181, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -5.52,': None,
'21181, 21181, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -5.999,': None,
'21550, 21181, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -5.52,': None,
'82132, 83661, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], 0.097,': None,
'88347, 93216, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.069,': None,
'79281, 58182, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -1.257,': None,
'21157, 19552, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -0.65,': None,
'88346, 94142, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.117,': None,
'57884, 59931, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -1.473,': None,
'01494, 00530, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -1.042,': None,
'99785, 92755, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -1.57,': None,
'17063, 24247, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.751,': None,
'16948, 17281, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -1.046,': None,
'92124, 94061, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -1.644,': None,
'19023, 18922, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.105,': None,
'13801, 13590, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2,': None,
'29390, 29447, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], 0.671,': None,
'60826, 92119, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -3.942,': None,
'09859, 10707, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.909,': None,
'79325, 81431, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -1.06,': None,
'79325, 81446, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -0.89,': None,
'93410, 04047, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.372,': None,
'17065, 16908, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.469,': None,
'04967, 04965, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.736,': None,
'83651, 93159, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], 2.342,': None,
'18850, 22119, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.276,': None,
'80423, 82835, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -0.974,': None,
'66590, 8639, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], 0.39,': None,
'82240, 83794, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -1.385,': None,
'80424, 82835, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -1.87,': None,
'94159, 97578, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -0.522,': None,
'94159, 99904, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -0.395,': None,
'00042, 01800, [1CH4], [1OH], [1:C3], [1:O3], [1:C3], [1:O3], -2.647,': None,
'00042, 01899, [1CH4], [1SH], [1:C3], [1:S3], [1:C3], [1:S3], -2.647,': None
       }

        # A few lines of CHEMBL data taken from test_dataset_goldeninput_pairs - ID's are chembl compound id's
        # the first twos lines represent the same pair but in both directions
        # a rest of the lines are a set of further unique CHEMBL ids taken from test_dataset_goldeninput_pairs but
        # fragments changed to same as the above to make more of the same pair change, should get 5 the same 1 different

        self.test_dataset_goldeninput_diffheader_02 = \
              'CUT,MOLID_L,MOLID_R,CONTEXT,FRAG_L,FRAG_R,ATTCHPT_CTX_L,ATTCHPT_FRAG_L,ATTCHPT_CTX_R,ATTCHPT_FRAG_R,ACT_A_DIFF,MOL_CLP_DIFF,MOL_MWT_DIFF'
      
        self.test_dataset_goldeninput_diffdata_02 = {
'single,205156,205260,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1cc[1cH]cc1,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],10.0,0.72,60.9': None,
'single,205260,205156,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1c[1cH]ccc1,Cc1cc[1cH]cc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],-10.0,-0.72,-60.9': None,
'single,383523,205156,[nH]1[n]c(c(c1)c1cc([n]cc1)c1cc[1cH]cc1)c1[n]cccc1,Cc1cc[1cH]cc1,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],11.0,0.72,60.9': None,
'single,382466,205260,[nH]1[n]c(c(c1)c1cc([n]cc1)c1c[1cH]ccc1)c1[n]cccc1,Cc1cc[1cH]cc1,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],8.0,0.72,60.9': None,
'single,382466,383433,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1cc[1cH]cc1,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],9.0,0.72,60.9': None,
'single,383433,426852,[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1cc[1cH]cc1,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],[1:CAR],[1:CAR],9.0,0.72,60.9': None
        }

        #
        # A matched pair taken from above data 205156, 205260 = CHEMBL205156 and CHEMBL205260)
        # Further lines have fake meaningless identifiers
        # self.test_dataset_goldeninput_diffheader_06
        #
        self.test_dataset_goldeninput_diffheader_03 = \
  'MOLID_L, MOLID_R, CONTEXT, FRAG_L, FRAG_R, ATTCHPT_CTX_L, ATTCHPT_FRAG_L, ATTCHPT_CTX_R, ATTCHPT_FRAG_R, DIFF_DIFF, MOL_CLP_DIFF, MOL_MWT_DIFF, MOL_L_CLP, MOL_R_CLP'

        self.test_dataset_goldeninput_diffdata_03 = {
  '205156, 205260, [nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1, Cc1cc[1cH]cc1, Cc1c[1cH]ccc1, [1:CAR], [1:CAR], [1:CAR], [1:CAR], -0.0628419184837, -0.441, -55.1, 0.000, 0.000': None,
  '9000001, 9000002, [nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1, Cc1cc[1cH]cc1, Cc1c[1cH]ccc1, [1:CAR], [1:CAR], [1:CAR], [1:CAR], -0.485163651843, -0.441, -55.1, 0.000, 0.000': None,
  '9000003, 9000004, [nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1, Cc1cc[1cH]cc1, Cc1c[1cH]ccc1, [1:CAR], [1:CAR], [1:CAR], [1:CAR], -0.183073670609, -0.441, -55.1, 0.000, 0.000': None,
  '9000005, 9000006, [nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1, Cc1cc[1cH]cc1, Cc1c[1cH]ccc1, [1:CAR], [1:CAR], [1:CAR], [1:CAR], -0.0826401755214, -0.68, -55.1, 0.000, 0.000': None,
  '9000007, 9000008, [nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1, Cc1cc[1cH]cc1, Cc1c[1cH]ccc1, [1:CAR], [1:CAR], [1:CAR], [1:CAR], -0.721352311615, -0.412, -55.1, 0.000, 0.000': None
        }

        #################################
        #
        # Golden results data
        #
        #################################

        # example aggregated results file
        self.test_dataset_goldenoutput_01 = {'FRAG_L,FRAG_R,CONTEXT_mean,CONTEXT_std,CONTEXT_count': None,
                                             '[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1cc[1cH]cc1,205260,NaN,1': None,
                                             '[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1c[1cH]ccc1,205156,NaN,1': None,
                                             '[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,Cc1[1cH]cccc1,205156,NaN,1': None,
                                             '[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1cc[1cH]cc1)(F)F,205156,NaN,1': None,
                                             '[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1c[1cH]ccc1)(F)F,205156,NaN,1': None,
                                             '[nH]1[n]c(c(c1)c1c[1cH][n]cc1)c1[n]cccc1,FC(c1[1cH]cccc1)(F)F,205156,NaN,1': None,
                                             '[nH]1[n]c(c(c1)c1cc([n]cc1)c1cc[1cH]cc1)c1[n]cccc1,[1CH4],383523,NaN,1': None,
                                             '[nH]1[n]c(c(c1)c1cc([n]cc1)c1cc[1cH]cc1)c1[n]cccc1,F[1CH](F)F,205156,NaN,1': None,
                                             '[nH]1[n]c(c(c1)c1cc([n]cc1)c1c[1cH]ccc1)c1[n]cccc1,[1CH4],382466,NaN,1': None,
                                             '[nH]1[n]c(c(c1)c1cc([n]cc1)c1c[1cH]ccc1)c1[n]cccc1,F[1CH](F)F,205260,NaN,1': None,
                                             '[nH]1[n]c(c(c1)c1cc(c2[1cH]cccc2)[n]cc1)c1[n]cccc1,[1CH4],426852,NaN,1': None,
                                             '[nH]1[n]c(c(c1)c1cc(c2[1cH]cccc2)[n]cc1)c1[n]cccc1,F[1CH](F)F,383433,NaN,1': None}

        self.test_dataset_goldenoutput_02 = {
            'FRAG_L,FRAG_R,ATTCHPT_CTX_L,ATTCHPT_CTX_R,PIC50_mean,PIC50_std,PIC50_count,CHEMBL_TARGET_ID_mean,CHEMBL_TARGET_ID_std,CHEMBL_TARGET_ID_count': None,
            '[1CH3]C,[1CH4],[1:C3],[1:C3],-0.974,0.076,2,0.000,0.000,2': None,
            '[1CH3]C,[1H],[1:C3],[1:C3],-2.339,0.096,2,0.000,0.000,2': None,
            '[1CH3]C,[1cH]1ccccc1,[1:C3],[1:C3],-2.545,0.146,2,0.000,0.000,2': None,
            '[1CH4],[1CH3]C,[1:C3],[1:C3],0.974,0.076,2,0.000,0.000,2': None,
            '[1CH4],[1H],[1:C3],[1:C3],-1.170,0.250,4,0.000,0.000,4': None,
            '[1CH4],[1cH]1ccccc1,[1:C3],[1:C3],-1.571,0.221,2,0.000,0.000,2': None,
            '[1H],[1CH3]C,[1:C3],[1:C3],2.339,0.096,2,0.000,0.000,2': None,
            '[1H],[1CH4],[1:C3],[1:C3],1.170,0.250,4,0.000,0.000,4': None,
            '[1H],[1cH]1ccccc1,[1:C3],[1:C3],-0.206,0.049,2,0.000,0.000,2': None,
            '[1cH]1ccccc1,[1CH3]C,[1:C3],[1:C3],2.545,0.146,2,0.000,0.000,2': None,
            '[1cH]1ccccc1,[1CH4],[1:C3],[1:C3],1.571,0.221,2,0.000,0.000,2': None,
            '[1cH]1ccccc1,[1H],[1:C3],[1:C3],0.206,0.049,2,0.000,0.000,2': None,
            'Cc1cc[1cH]cc1,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],-0.834,NaN,1,0.000,NaN,1': None,
            'Cc1cc[1cH]cc1,Cc1[1cH]cccc1,[1:CAR],[1:CAR],-2.153,NaN,1,0.000,NaN,1': None,
            'Cc1cc[1cH]cc1,FC(c1cc[1cH]cc1)(F)F,[1:CAR],[1:CAR],0.067,NaN,1,0.000,NaN,1': None,
            'Cc1cc[1cH]cc1,FC(c1c[1cH]ccc1)(F)F,[1:CAR],[1:CAR],-0.933,NaN,1,0.000,NaN,1': None,
            'Cc1cc[1cH]cc1,FC(c1[1cH]cccc1)(F)F,[1:CAR],[1:CAR],-2.523,NaN,1,0.000,NaN,1': None,
            'Cc1c[1cH]ccc1,Cc1cc[1cH]cc1,[1:CAR],[1:CAR],0.834,NaN,1,0.000,NaN,1': None,
            'Cc1c[1cH]ccc1,Cc1[1cH]cccc1,[1:CAR],[1:CAR],-1.319,NaN,1,0.000,NaN,1': None,
            'Cc1c[1cH]ccc1,FC(c1cc[1cH]cc1)(F)F,[1:CAR],[1:CAR],0.901,NaN,1,0.000,NaN,1': None,
            'Cc1c[1cH]ccc1,FC(c1c[1cH]ccc1)(F)F,[1:CAR],[1:CAR],-0.099,NaN,1,0.000,NaN,1': None,
            'Cc1c[1cH]ccc1,FC(c1[1cH]cccc1)(F)F,[1:CAR],[1:CAR],-1.689,NaN,1,0.000,NaN,1': None,
            'Cc1[1cH]cccc1,Cc1cc[1cH]cc1,[1:CAR],[1:CAR],2.153,NaN,1,0.000,NaN,1': None,
            'Cc1[1cH]cccc1,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],1.319,NaN,1,0.000,NaN,1': None,
            'Cc1[1cH]cccc1,FC(c1cc[1cH]cc1)(F)F,[1:CAR],[1:CAR],2.220,NaN,1,0.000,NaN,1': None,
            'Cc1[1cH]cccc1,FC(c1c[1cH]ccc1)(F)F,[1:CAR],[1:CAR],1.220,NaN,1,0.000,NaN,1': None,
            'Cc1[1cH]cccc1,FC(c1[1cH]cccc1)(F)F,[1:CAR],[1:CAR],-0.370,NaN,1,0.000,NaN,1': None,
            'FC(c1cc[1cH]cc1)(F)F,Cc1cc[1cH]cc1,[1:CAR],[1:CAR],-0.067,NaN,1,0.000,NaN,1': None,
            'FC(c1cc[1cH]cc1)(F)F,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],-0.901,NaN,1,0.000,NaN,1': None,
            'FC(c1cc[1cH]cc1)(F)F,Cc1[1cH]cccc1,[1:CAR],[1:CAR],-2.220,NaN,1,0.000,NaN,1': None,
            'FC(c1cc[1cH]cc1)(F)F,FC(c1c[1cH]ccc1)(F)F,[1:CAR],[1:CAR],-1.000,NaN,1,0.000,NaN,1': None,
            'FC(c1cc[1cH]cc1)(F)F,FC(c1[1cH]cccc1)(F)F,[1:CAR],[1:CAR],-2.590,NaN,1,0.000,NaN,1': None,
            'FC(c1c[1cH]ccc1)(F)F,Cc1cc[1cH]cc1,[1:CAR],[1:CAR],0.933,NaN,1,0.000,NaN,1': None,
            'FC(c1c[1cH]ccc1)(F)F,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],0.099,NaN,1,0.000,NaN,1': None,
            'FC(c1c[1cH]ccc1)(F)F,Cc1[1cH]cccc1,[1:CAR],[1:CAR],-1.220,NaN,1,0.000,NaN,1': None,
            'FC(c1c[1cH]ccc1)(F)F,FC(c1cc[1cH]cc1)(F)F,[1:CAR],[1:CAR],1.000,NaN,1,0.000,NaN,1': None,
            'FC(c1c[1cH]ccc1)(F)F,FC(c1[1cH]cccc1)(F)F,[1:CAR],[1:CAR],-1.590,NaN,1,0.000,NaN,1': None,
            'FC(c1[1cH]cccc1)(F)F,Cc1cc[1cH]cc1,[1:CAR],[1:CAR],2.523,NaN,1,0.000,NaN,1': None,
            'FC(c1[1cH]cccc1)(F)F,Cc1c[1cH]ccc1,[1:CAR],[1:CAR],1.689,NaN,1,0.000,NaN,1': None,
            'FC(c1[1cH]cccc1)(F)F,Cc1[1cH]cccc1,[1:CAR],[1:CAR],0.370,NaN,1,0.000,NaN,1': None,
            'FC(c1[1cH]cccc1)(F)F,FC(c1cc[1cH]cc1)(F)F,[1:CAR],[1:CAR],2.590,NaN,1,0.000,NaN,1': None,
            'FC(c1[1cH]cccc1)(F)F,FC(c1c[1cH]ccc1)(F)F,[1:CAR],[1:CAR],1.590,NaN,1,0.000,NaN,1': None,
            '[1CH4],F[1CH](F)F,[1:CAR],[1:CAR],-0.134,0.221,3,0.000,0.000,3': None,
            'F[1CH](F)F,[1CH4],[1:CAR],[1:CAR],0.134,0.221,3,0.000,0.000,3': None}

        # Note that for this output set the N value is 1 (not 4) due to the duplicate removal we applied
        # so for fragment change [2CH3][1CH3] --> [12CH4] we only count it once! not 4 times for the same 
        # id/id pair
        self.test_dataset_goldenoutput_03 = {'FRAG_L,FRAG_R,PIC50_DIFF_mean,PIC50_DIFF_std,PIC50_DIFF_count': None,
                                             '[2CH3][1CH3],[12CH4],-0.100,NaN,1': None,
                                             '[1CH3][2CH3],[12CH4],-0.100,NaN,1': None,
                                             '[12CH4],[2CH3][1CH3],0.100,NaN,1': None,
                                             '[12CH4],[1CH3][2CH3],0.100,NaN,1': None,
                                             '[2CH3]C[1CH3],[2CH3][1CH3],-0.100,NaN,1': None,
                                             '[1CH3]C[2CH3],[1CH3][2CH3],-0.100,NaN,1': None,
                                             '[2CH3][1CH3],[2CH3]C[1CH3],0.100,NaN,1': None,
                                             '[1CH3][2CH3],[1CH3]C[2CH3],0.100,NaN,1': None}

        # same as self.test_dataset_goldenoutput_03 but the counts 
        # are higher due to id dupes for given fragL, fragR change (no filtering)
        self.test_dataset_goldenoutput_03b = {'FRAG_L,FRAG_R,PIC50_DIFF_mean,PIC50_DIFF_std,PIC50_DIFF_count': None,
                                              '[2CH3][1CH3],[12CH4],-0.100,0.000,4': None,
                                              '[1CH3][2CH3],[12CH4],-0.100,0.000,4': None,
                                              '[12CH4],[2CH3][1CH3],0.100,0.000,4': None,
                                              '[12CH4],[1CH3][2CH3],0.100,0.000,4': None,
                                              '[2CH3]C[1CH3],[2CH3][1CH3],-0.100,NaN,1': None,
                                              '[1CH3]C[2CH3],[1CH3][2CH3],-0.100,NaN,1': None,
                                              '[2CH3][1CH3],[2CH3]C[1CH3],0.100,NaN,1': None,
                                              '[1CH3][2CH3],[1CH3]C[2CH3],0.100,NaN,1': None}
        
        # diff function test output
        self.test_dataset_goldenoutput_04 = {
         'Mean, SD, DIFF60, upp90, low90': None,
         '51.0, -1.60, -56.4, -51.5, -58.5': None
        }
        
        self.test_dataframe_goldenoutput_04 = {
          'FRAG_L,FRAG_R,DIFF_DIFF_count,DIFF_DIFF_n_pos_diff,DIFF_DIFF_n_neg_diff,DIFF_DIFF_diff60_low,DIFF_DIFF_diff60,DIFF_DIFF_diff60_upp': None,
          '[1CH4],[1SH],1,0.000,1.000,NaN,-59.663,NaN': None,
          '[1CH4],[1OH],50,8.000,42.000,-58.6488,-56.7383,-52.336': None
        }

        self.test_dataframe_goldenoutput_05 = {
            'FRAG_L,FRAG_R,ACT_A_DIFF_diff60_low,ACT_A_DIFF_diff60,ACT_A_DIFF_diff60_upp,ACT_A_DIFF_count,ACT_A_DIFF_n_pos_diff,ACT_A_DIFF_n_neg_diff,MOL_CLP_DIFF_mean,MOL_CLP_DIFF_std,MOL_MWT_DIFF_first': None,
            'Cc1cc[1cH]cc1,Cc1c[1cH]ccc1,40.0,40.0,40.0,5,5.000,0.000,0.720,0.000,60.900': None,
            'Cc1c[1cH]ccc1,Cc1cc[1cH]cc1,NaN,-60.0,NaN,1,0.000,1.000,-0.720,NaN,-60.900': None
        }

        self.test_dataframe_goldenoutput_06 = {
            'FRAG_L,FRAG_R,DIFF_DIFF_mean_diff_invlog_low,DIFF_DIFF_mean_diff_invlog,DIFF_DIFF_mean_diff_invlog_upp,DIFF_DIFF_count,DIFF_DIFF_n_pos_diff,DIFF_DIFF_n_neg_diff,MOL_CLP_DIFF_mean,MOL_CLP_DIFF_std,MOL_MWT_DIFF_first': None,
            'Cc1cc[1cH]cc1,Cc1c[1cH]ccc1,0.2628,0.4932,0.9253,5,0.000,5.000,-0.483,0.111,-55.100': None
        }

        self.test_dataframe_goldenoutput_07 = {
            'FRAG_L,FRAG_R,DIFF_DIFF_mean_diff_low,DIFF_DIFF_mean_diff,DIFF_DIFF_mean_diff_upp,DIFF_DIFF_count,DIFF_DIFF_n_pos_diff,DIFF_DIFF_n_neg_diff,MOL_CLP_DIFF_mean,MOL_CLP_DIFF_std,MOL_MWT_DIFF_first': None,
            'Cc1cc[1cH]cc1,Cc1c[1cH]ccc1,-0.5803,-0.307,-0.0337,5,0.000,5.000,-0.483,0.111,-55.100': None
        }

        self.test_dataframe_goldenoutput_08 = {
            'FRAG_L,FRAG_R,DIFF_DIFF_mean_diff_low,DIFF_DIFF_mean_diff,DIFF_DIFF_mean_diff_upp,DIFF_DIFF_count,DIFF_DIFF_n_pos_diff,DIFF_DIFF_n_neg_diff,MOL_CLP_DIFF_mean,MOL_CLP_DIFF_std,MOL_MWT_DIFF_first,MOL_L_CLP_std,QUALITY': None,
            'Cc1cc[1cH]cc1,Cc1c[1cH]ccc1,-0.5803,-0.307,-0.0337,5,0.000,5.000,-0.483,0.111,-55.100,0.000,0': None
        }

        ##################################
        #
        # write test data to temp files
        #
        ##################################
        
        # csv file
        self.temp_file_input_csv.write(', '.join(self.test_dataset_goldeninput_csv_headers)+"\n")
        for data in list(self.test_dataset_goldeninput_csv_data.keys()):
            self.temp_file_input_csv.write(data+"\n")
        self.temp_file_input_csv.close()

        # std pairs
        self.temp_file_input_pairs.write(self.test_dataset_goldeninput_header+"\n")

        for line in self.test_dataset_goldeninput_pairs.items():
            #print(line)
            self.temp_file_input_pairs.write(line[0]+"\n")
        self.temp_file_input_pairs.close()

        # pairs after de-dupe
        self.temp_file_input_pairdupes.write(self.test_dataset_goldeninput_pairdupes_header+"\n")
        for line in self.test_dataset_goldeninput_pairdupes.items():
            self.temp_file_input_pairdupes.write(line[0]+"\n")

        self.temp_file_input_pairdupes.close()

        # diff function test
        self.temp_file_input_diffdata.write(self.test_dataset_goldeninput_diffheader+"\n")
        for line in self.test_dataset_goldeninput_diffdata.items():
            self.temp_file_input_diffdata.write(line[0]+"\n")

        self.temp_file_input_diffdata.close()

        # diff function test with props
        self.temp_file_input_diffdata_02.write(self.test_dataset_goldeninput_diffheader_02+"\n")
        for line in self.test_dataset_goldeninput_diffdata_02.items():
            self.temp_file_input_diffdata_02.write(line[0]+"\n")

        self.temp_file_input_diffdata_02.close()
        
        # container for results data
        self.test_dataset_testresults = {}

    def tearDown(self):        
        """Tear down object for clean reuse in further tests"""

        self.test_dataset_testresults.clear()
        
        os.remove(self.temp_file_input_pairs.name)
        os.remove(self.temp_file_input_diffdata.name)
        os.remove(self.temp_file_input_pairdupes.name)

    def test_diff60_df(self):
        """Functional test to confirm the diff function is working on a dataframe"""
        
        # group, aggregate, convert object to df, sort index
        grouped = self.test_dataframe.groupby(['A'])
        grouped = grouped['C'].agg([np.mean, len, mmps.diffn_agg(60)])
        
        # this is important for unit testing as order or columns can change and needs to be made consistent
        grouped = pd.DataFrame(grouped).reset_index()
        grouped.sort_index(axis=1, inplace=True)

        #print(self.test_dataframe_result)
        pd.util.testing.assert_frame_equal(grouped, self.test_dataframe_result)

    def test_pairsdataobj_to_pd(self):
        """Test method to create a dataframe from the base object iterator
        This test and its associated method are the only ones that use the methods and 
        other objects associated with the class we inherit from (MMPDataObjectsClass)
        """

        # full build of pairs and data objects from csv
        self.test_mmp_pairs_object.csv_sniffer(self.temp_file_input_csv.name, 'SMILES', 'CHEMBL_ID')
        self.test_mmp_pairs_object.csv_to_data_objects(self.temp_file_input_csv.name, 'SMILES', 'CHEMBL_ID')
        tmp_dicer_file = self.test_mmp_pairs_object.write_mol_smi_dict_tofile()
        self.test_mmp_pairs_object.build_from_dicer(tmp_dicer_file, 'BOTH', 'NONE')
        # then convert final pairs to pandas dataframe
        self.test_mmp_pairs_object.pairsdataobj_to_pd('BOTH', 'PIC50')

        columns = list(self.test_mmp_pairs_object.pairs_table.columns.values)
        shape = self.test_mmp_pairs_object.pairs_table.shape
        num_nulls = self.test_mmp_pairs_object.pairs_table.isnull().values.sum()
        # don't want whole table in unittest so just check structure
        self.assertEqual(columns,
                         ['CUT', 'MOLID_L', 'MOLID_R', 'CONTEXT', 'FRAG_L', 'FRAG_R', 'ATTCHPT_CTX_L', 'ATTCHPT_FRAG_L',
                          'ATTCHPT_CTX_R', 'ATTCHPT_FRAG_R', 'PIC50', 'CHEMBL_TARGET_ID'])
        self.assertEqual(shape[0], 280)
        self.assertEqual(shape[1], 12)
        self.assertEqual(num_nulls, 0)

    def test_pd_aggregation_frag_filter_dupes(self):
        """Test method to filter certain dupes from a dataframe"""

        self.test_mmp_pairs_object.pd_read_csv(self.temp_file_input_pairdupes.name, self.test_columns_dict)
        self.test_mmp_pairs_object.pd_aggregate_pairs_to_csv(self.temp_file_output_nodupes.name, 'frag')

        test_results_filehandle = open(self.temp_file_output_nodupes.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None

        #print(self.test_dataset_testresults)
        self.assertEqual(self.test_dataset_testresults, self.test_dataset_goldenoutput_03)

    def test_pd_aggregation_frag_filter_dupes_false(self):
        """Test method to filter certain dupes from a dataframe"""

        self.test_mmp_pairs_object.pd_read_csv(self.temp_file_input_pairdupes.name, self.test_columns_dict)
        self.test_mmp_pairs_object.pd_aggregate_pairs_to_csv(self.temp_file_output_nodupes.name, 'frag',
                                                             remove_id_dupes=False)
        
        test_results_filehandle = open(self.temp_file_output_nodupes.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None

        #print(self.test_dataset_testresults)
        self.assertEqual(self.test_dataset_testresults, self.test_dataset_goldenoutput_03b)

    def test_pd_aggregation_frag(self):
        """Build pandas df from csv, aggregate then save pairs"""

        # do everything, read, group by frag (not frag attach), write
        self.test_mmp_pairs_object.pd_read_csv(self.temp_file_input_pairs.name, self.test_columns_dict)
        self.test_mmp_pairs_object.pd_aggregate_pairs_to_csv(self.temp_file_output.name, 'frag')

        test_results_filehandle = open(self.temp_file_output.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None

        #print(self.test_dataset_testresults)
        self.assertEqual(self.test_dataset_testresults, self.test_dataset_goldenoutput_01)

    def test_pd_aggregation_attach(self):
        """Build pandas df from mmp data object iterator, aggregate then save pairs"""

        # self.test_dataset_goldeninput_csv_data
        # full build of pairs and data objects from csv
        self.test_mmp_pairs_object.csv_sniffer(self.temp_file_input_csv.name, 'SMILES', 'CHEMBL_ID')
        self.test_mmp_pairs_object.csv_to_data_objects(self.temp_file_input_csv.name, 'SMILES', 'CHEMBL_ID')
        tmp_dicer_file = self.test_mmp_pairs_object.write_mol_smi_dict_tofile()
        self.test_mmp_pairs_object.build_from_dicer(tmp_dicer_file, 'SINGLE', 'NONE')
        # then convert final pairs to pandas dataframe
        self.test_mmp_pairs_object.pairsdataobj_to_pd('SINGLE', 'PIC50')

        self.test_mmp_pairs_object.pd_aggregate_pairs_to_csv(self.temp_file_output.name, 'attach')

        test_results_filehandle = open(self.temp_file_output.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None

        #print(self.test_dataset_testresults) # use this pprint statement to regenerate the golden data
        self.assertEqual(self.test_dataset_testresults, self.test_dataset_goldenoutput_02)

    def test_pd_aggregation_diff(self):
        """Test the diff function"""

        # do everything, read, group by frag (not frag attach), write
        self.test_mmp_pairs_object.pd_read_csv(self.temp_file_input_diffdata.name)
        self.test_mmp_pairs_object.pd_aggregate_pairs_to_csv(self.temp_file_output.name,
                                                             'frag',
                                                             agg_method='diff60',
                                                             act_col='DIFF')
        # self.test_mmp_pairs_object.grouped = pd.DataFrame(self.test_mmp_pairs_object.grouped).reset_index()
        # self.test_mmp_pairs_object.grouped.sort_index(axis=1, inplace=True)
        # self.test_mmp_pairs_object.pd_save_aggregate_pairs(self.temp_file_output.name)

        test_results_filehandle = open(self.temp_file_output.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None

        #print(self.test_dataset_testresults)
        self.assertEqual(self.test_dataframe_goldenoutput_04, self.test_dataset_testresults)

    def test_pd_aggregation_diff_withprops(self):
        """Test the diff function with prop data in same dataframe"""

        # do everything, read, group by frag (not frag attach), write
        self.test_mmp_pairs_object.pd_read_csv(self.temp_file_input_diffdata_02.name)
        self.test_mmp_pairs_object.pd_aggregate_pairs_to_csv(self.temp_file_output.name,
                                                             'frag',
                                                             agg_method='diff60',
                                                             prop_data=True,
                                                             act_col='ACT_A')
        # now read results back in for unittest compare
        test_results_filehandle = open(self.temp_file_output.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None

        print(">>05>>", self.test_dataset_testresults)
        self.assertEqual(self.test_dataframe_goldenoutput_05, self.test_dataset_testresults)

    def test_pd_aggregation_meandiff_invlog_withprops(self):

        self.temp_file_input_meandiff_withprops = tempfile.NamedTemporaryFile(delete=False,
                                                                              encoding='utf-8',
                                                                              mode='wt')
        
        # diff function test with props
        self.temp_file_input_meandiff_withprops.write(self.test_dataset_goldeninput_diffheader_03+"\n")
        for line in self.test_dataset_goldeninput_diffdata_03.items():
            self.temp_file_input_meandiff_withprops.write(line[0]+"\n")
        self.temp_file_input_meandiff_withprops.close()
                
        self.test_mmp_pairs_object.pd_read_csv(self.temp_file_input_meandiff_withprops.name)
        self.test_mmp_pairs_object.pd_aggregate_pairs_to_csv(self.temp_file_output.name,
                                                             'frag',
                                                             agg_method='mean_diff_invlog',
                                                             prop_data=True,
                                                             act_col='DIFF')

        test_results_filehandle = open(self.temp_file_output.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None

        os.remove(self.temp_file_input_meandiff_withprops.name)

        #print(">>06>>", self.test_dataset_testresults)
        self.assertEqual(self.test_dataframe_goldenoutput_06, self.test_dataset_testresults)

    def test_pd_aggregation_meandiff_withprops(self):
        """ """
        self.temp_file_input_meandiff_withprops = tempfile.NamedTemporaryFile(delete=False,
                                                                              encoding='utf-8',
                                                                              mode='wt')

        # diff function test with props
        self.temp_file_input_meandiff_withprops.write(self.test_dataset_goldeninput_diffheader_03+"\n")
        for line in self.test_dataset_goldeninput_diffdata_03.items():
            self.temp_file_input_meandiff_withprops.write(line[0]+"\n")
        self.temp_file_input_meandiff_withprops.close()

        self.test_mmp_pairs_object.pd_read_csv(self.temp_file_input_meandiff_withprops.name)
        self.test_mmp_pairs_object.pd_aggregate_pairs_to_csv(self.temp_file_output.name,
                                                             'frag',
                                                             agg_method='mean_diff',
                                                             prop_data=True,
                                                             act_col='DIFF')

        test_results_filehandle = open(self.temp_file_output.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None

        os.remove(self.temp_file_input_meandiff_withprops.name)

        #print(">>07>>", self.test_dataset_testresults)
        self.assertEqual(self.test_dataframe_goldenoutput_07, self.test_dataset_testresults)

    def test_pd_aggregation_meandiff_withprops_withqmetric(self):
        """Test of the Quality metric, output data last column should be 0,1 or 2 representing quality metric"""
        self.temp_file_input_meandiff_withprops = tempfile.NamedTemporaryFile(delete=False,
                                                                              encoding='utf-8',
                                                                              mode='wt')

        # diff function test with props
        self.temp_file_input_meandiff_withprops.write(self.test_dataset_goldeninput_diffheader_03+"\n")
        for line in self.test_dataset_goldeninput_diffdata_03.items():
            self.temp_file_input_meandiff_withprops.write(line[0]+"\n")
        self.temp_file_input_meandiff_withprops.close()


        self.test_mmp_pairs_object.pd_read_csv(self.temp_file_input_meandiff_withprops.name)
        self.test_mmp_pairs_object.pd_aggregate_pairs_to_csv(self.temp_file_output.name,
                                                             'frag',
                                                             agg_method='mean_diff',
                                                             prop_data=True,
                                                             act_col='DIFF',
                                                             add_qmetric=True)

        test_results_filehandle = open(self.temp_file_output.name, 'r')
        for line in test_results_filehandle:
            line = line.rstrip('\r')
            line = line.rstrip('\n')
            self.test_dataset_testresults[line] = None

        os.remove(self.temp_file_input_meandiff_withprops.name)

        print(self.test_dataset_testresults)
        self.assertEqual(self.test_dataframe_goldenoutput_08, self.test_dataset_testresults)


if __name__ == '__main__':
    unittest.main()
