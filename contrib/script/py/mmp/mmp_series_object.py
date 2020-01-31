###################################################################
""" Summary: Class and Methods to get and interogate matched series

About: An object that reads a smiles file into memory and generates
matched series (OBoyle, Bostron, Sayle, Gill JMC 2014) of a defined
length via an in-house extension of the fragment indexing method
(Hussain and Rea JCIM 2010). Includes various matched series scoring
methods as annotated. (JAL)

TODO: print Molecular Context from result series.  Currently the context of the query series is printed.
TODO: print Series source as filename of query series or blank for internal SARA transfer

"""
###################################################################
import logging
import sys
import os
import glob
from copy import deepcopy

# specific
import pandas as pd
from operator import itemgetter
# from pandas.util.testing import assert_frame_equal

#
from itertools import groupby, combinations  # product, permutations #
from math import factorial
from scipy.stats import binom_test, pearsonr, spearmanr, linregress, skew
from scipy.spatial.distance import cityblock
import numpy as np
# from numpy import NaN

# local imports
from mmp.mmp_data_objects import MMPDataObjectClass
from mmp.mmp_math_functions import inv_cantor
import mmp.mmp_enum_mols_from_pairs as enum_mols

# things needed for unit testing only
import unittest
import tempfile
#from tempfile import NamedTemporaryFile, mkdtemp

# https://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
pd.options.mode.chained_assignment = None


class MMPSeriesObjectClass(MMPDataObjectClass):

    """Class implements objects and methods for MMP series manipulation"""

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

        # variable use to track method used to create series
        self.data_column = None
        self.data_column_position = None
        self.min_series_length = 3
        self.return_series_max_len = None
        self.threshold = 0.50001
        # series number incremented by _iterators
        self.series_id = 0
        # filtering out series:
        self.series_filtered_out = 0
        self.max_skew = 3
        self.min_pAct = 0.5
        #
        self.running_fullscan = False

        # base pandas df
        self.series_df = None
        self.ref_series_df = None

        # temp stuff
        self.key_smi_tempfile = None
        self.series_comparison_df = None
        self.enumerated_products_smi = None
        # this becomes series context => max activity
        self.max_series_potency = {}

        # set of series files to inject
        self.series_from_file = False
        self.series_base_dir = None
        self.series_file_list = []

    def clean_out_data_seriesobj(self):
        """Implements method from MMPDataObjectClass to clear out all data structures plus
        clears out a few more related to this extended class"""

        self.clean_out_data()

        self.data_column = None
        self.data_column_position = None
        self.min_series_length = 3
        self.return_series_max_len = None
        self.threshold = 0.50001
        self.series_id = 0
        #
        self.series_filtered_out = 0
        self.max_skew = 3
        self.min_pAct = 0.5

        self.series_df = None
        self.ref_series_df = None

        self.key_smi_tempfile = None
        self.series_comparison_df = None
        self.enumerated_products_smi = None
        self.max_series_potency = {}

        self.series_from_file = False
        self.series_base_dir = None
        self.series_file_list = []

    ############################################################################################
    #
    # Setup and series retrieval methods
    #
    ############################################################################################

    def setup_mmp_data_for_mms(self,
                               input_file,
                               smi_col,
                               id_col,
                               data_column,
                               min_series_len,
                               threshold,
                               cut_type='BOTH'):
        """
        Method sets object level vars for MMS generation as well as running a number setup methods from mmp
        object to generate and store pairs ready for MMS generation.  Can take some time to do this.
        :param input_file: CSV of smiles, id and activity data
        :param smi_col: column name of smiles column
        :param id_col: column name of id column
        :param data_column: column name of activity data column
        :param min_series_len: minimum length of a series must be >3
        :param threshold: maxff dicer setting (% molecule to be retained as context)
        :param cut_type: single double or both cuts?
        :return: None
        """
        cut_type_options = ['SINGLE', 'DOUBLE', 'BOTH']
        if cut_type not in cut_type_options:
            raise Exception('Errro, invalid cut type: %s' % cut_type)
        self.logger.info('Setting up MMS generation via MMP data creation')
        # parse csv
        self.csv_sniffer(input_file, smi_col, id_col)
        self.csv_to_data_objects(input_file, smi_col, id_col, std_smi=True)
        # get pairs
        tmp_dicer_file = self.write_mol_smi_dict_tofile()
        ##############################################################################
        # This line check for odd stuff in smiles like chirality and fails if present
        ##############################################################################
        self.scan_input_smiles(tmp_dicer_file, fail_chirals=True)
        # build series
        self.build_from_dicer(tmp_dicer_file, cut_type, 'NONE', threshold=threshold)
        #
        self.min_series_length = min_series_len
        self.data_column = data_column
        self.threshold = threshold
        #
        # Find the column position for the data_column in self.mol_data_dict
        try:
            self.data_column_position = self.headers_nosmi.index(data_column)
        except:
            raise Exception("Invalid data column specification, %s cannot be found in CSV headers" % data_column)

        self.logger.info('completed Setting up MMP data: %d sgl & %d dbl ctx stored for series gen' %
                         (len(self.single_pairs_dict), len(self.double_pairs_dict)))

    def _iterator_mmp_series_numeric(self,
                                     use_comparison_df=False,
                                     store_series_max_act=True,
                                     sgl_or_dbl='single',
                                     apply_pre_filter=False):
        """ Method to iterate over Single Cut Dictionary structure and yield seed matched series
        This does not exhaustively return all matched series
        """

        if self.min_series_length < 3:
            sys.exit("Invalid value for param min_series_length in method iterator_singlecut_mmp_series, must be > 3")

        #
        # can later change this to check which dict we are comparing, allows search between SMI sets not just
        # within SMI sets via use of comparison dict
        if sgl_or_dbl == 'single':
            if use_comparison_df:
                query_dict = self.single_pairs_comparison_dict
            else:
                query_dict = self.single_pairs_dict

            self.logger.info('Iterating over single cut pairs dictionary to get matched series of length N: %d'
                             % self.min_series_length)

        elif sgl_or_dbl == 'double':
            if use_comparison_df:
                query_dict = self.double_pairs_comparison_dict
            else:
                query_dict = self.double_pairs_dict

            self.logger.info('Iterating over double cut pairs dictionary to get matched series of length N: %d'
                             % self.min_series_length)

        else:
            self.logger.warn('invalid cut type: %s (try single or double)' % sgl_or_dbl)
            sys.exit('invalid cut type: %s (try single or double)' % sgl_or_dbl)

        #
        for ctx_id in query_dict:

            # series are already determined by the data structure, for example:
            #     dict_single[ctx1]  => { mol1_frag1_id => frag1, mol2_frag4_id => frag4, ... } is a series
            #     dict_single[ctx2]  => { mol1_frag2_id => frag2, mol2_frag5_id => frag5, ... } is a series
            #     dict_single[ctx3]  => { mol1_frag3_id => frag3, mol7_frag9_id => frag7} is only a pair
            #     dict_single[ctx4]  => { mol2_frag6_id => frag6 } is nothing!
            # series exist where the context is the same, but the fragments differ and there are more
            # than two of them i.e.: param series_len > 2
            # ...but complexity occurs where we have two different mol_id's with the same frag_id's but diff
            # activity. It also arises when two different mol_id's and frag_id's have identical activities

            series_length = len(query_dict[ctx_id])
            series_unsorted = []

            if series_length >= self.min_series_length:

                self.logger.debug('Found a series of valid length (%s)' % series_length)

                # iterate and store all fragments, create series
                for molid_fragid_uid_L in query_dict[ctx_id]:

                    molid_L, frag_id_L = inv_cantor(molid_fragid_uid_L)
                    # frag_L = self.refsmi_dict[frag_id_L]

                    # TODO: Code does not support situation where there are 1 or more text based columns!
                    # print self.mol_data_dict[molid_L]
                    # print self.mol_data_dict
                    # print self.headers
                    # print self.headers_nosmi
                    # print self.headers_numeric_position
                    try:
                        mol_act = self.mol_data_dict[molid_L][self.data_column_position]
                    except:
                        mol_act = None

                    series_unsorted.append([frag_id_L, molid_L, mol_act])

                #################################################################################
                # filter out series if:
                # (a) the set of unique activity values is less than self.min_series_length
                # or they have poor characteristics, Keefer & Chang MedChemComm 2017 - 10.1039/C7MD00465F
                # (b) range in activity <= 0.5
                # (c) skew <= 3
                # quit the loop, discarding the series
                #################################################################################
                act_data_set = set()
                act_data_arr = []
                for item in series_unsorted:
                    act_data_set.add(item[2])
                    act_data_arr.append(item[2])
                # print act_data_set
                # print len(act_data_set), self.min_series_length
                # (a)
                if len(act_data_set) < self.min_series_length:
                    continue

                # now sort it
                self.logger.debug('Ordering series')
                series_sorted = sorted(series_unsorted, key=lambda xx: xx[2])

                # TODO: check this as it seems too discriminatory!
                act_data_nparr = np.array(act_data_arr)
                # print ">>> range and skew"
                # print "data: ", act_data_nparr
                # print "range: ", (abs(act_data_nparr.max() - act_data_nparr.min()))
                # print "skew: ", (3 * (act_data_nparr.mean() - np.median(act_data_nparr)))/act_data_nparr.std()
                # print skew(act_data_nparr)
                if apply_pre_filter:

                    act_data_nparr = np.array(act_data_arr)

                    # remove if below pAct range
                    if (abs(act_data_nparr.max() - act_data_nparr.min())) <= self.min_pAct:
                        self.series_filtered_out += 1
                        continue

                    # skew, via pearsons second coefficient (alternate method avoiding mode)
                    # remove if series is skewed
                    # http://www.statisticshowto.com/pearson-mode-skewness/
                    # if (3 * (act_data_nparr.mean() - np.median(act_data_nparr)))/act_data_nparr.std() <= 3:
                    # easier using scipy:
                    # https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.stats.skew.html
                    if skew(act_data_nparr) >= self.max_skew:
                        self.series_filtered_out += 1
                        continue

                self.series_id += 1
                series_seq_id = 0
                previous_act = None

                if store_series_max_act:
                    self.max_series_potency[ctx_id] = series_sorted[-1][2]

                #####################################################################
                # A sorted series gets a series_id
                # Each fragment in the sorted series gets a series_seq_id
                # The series_seq_id will only be incremented if the fragment has a greater activity than the previous
                # activity.  Identical series_seq_id within a given series indicate identical activity fragments
                #####################################################################
                for [frag_id_L, molid_L, mol_act] in series_sorted:

                    if previous_act is None:
                        previous_act = mol_act
                        series_seq_id += 1

                    else:
                        if mol_act > previous_act:
                            series_seq_id += 1
                            previous_act = mol_act

                    # (ctx1_id, ctx2_id) = inv_cantor(ctx_id)
                    # print self.series_id, series_seq_id, self.refsmi_dict[ctx1_id], self.refsmi_dict[ctx2_id],
                    #   self.refsmi_dict[frag_id_L], molid_L, mol_act
                    yield self.series_id, series_seq_id, ctx_id, frag_id_L, molid_L, mol_act

        self.logger.info('Removed %s series due to skew and poor pAct range' % self.series_filtered_out)
        # print ('Removed %s series due to skew and poor pAct range' % self.series_filtered_out)

    def generate_store_mmp_series(self, sgl_or_dbl_or_both='single', apply_pre_filter=False):
        """
        Utilises MMS iterator to get single cut MMS's and store them
        :param sgl_or_dbl_or_both: single cut, double cut or both
        :param apply_pre_filter
        :return: nothing, updates self.series_df in place
        """

        series_data = []

        if sgl_or_dbl_or_both == 'single' or sgl_or_dbl_or_both == 'both':
            series_data_sgl = [(a, b, c, d, e, f) for a, b, c, d, e, f in
                               self._iterator_mmp_series_numeric(apply_pre_filter=apply_pre_filter)]
            series_data.extend(series_data_sgl)
            # print series_data

        if sgl_or_dbl_or_both == 'double' or sgl_or_dbl_or_both == 'both':
            series_data_dbl = [(a, b, c, d, e, f) for a, b, c, d, e, f in
                               self._iterator_mmp_series_numeric(sgl_or_dbl='double',
                                                                 apply_pre_filter=apply_pre_filter)]
            series_data.extend(series_data_dbl)
            # print series_data

        self.series_df = pd.DataFrame(series_data,
                                      columns=['SERIES_ID', 'SERIES_SEQ_ID',
                                               'CONTEXT_ID', 'FRAG_ID',
                                               'MOL_ID', 'ACTIVITY']
                                      )
        # self.series_df.set_index(['SERIES_ID', 'FRAG_ID'], inplace=True)

        # print('Parsed series CSV to dataframe of size %d, %d' % (self.series_df.shape[0], self.series_df.shape[1]))
        self.logger.info('Parsed series CSV to dataframe of size %d, %d' %
                         (self.series_df.shape[0], self.series_df.shape[1]))

    ############################################################################################
    #
    # Search Methods
    #
    ############################################################################################

    def return_series_matching_idlist(self, id_list, column_name, use_comparison_df=False, strict_ordering=False):
        """
        Search a df for all series that individually match every one of the frag ids in the query fragid list
        :param id_list: a list of fragment ID's to use as query
        :param column_name: the name of a column in self.series_df that will be searched for the ids listed
        :param use_comparison_df compare two dfs
        :param strict_ordering
        :return: A pandas df of series that match the query
        """
        if use_comparison_df:
            self.logger.info('Searching comparison_df for series matching id list %s of type %s' %
                             (str(id_list), column_name))
        else:
            self.logger.info('Searching series_df for series matching id list %s of type %s' %
                             (str(id_list), column_name))

        if len(id_list) < 2:
            raise Exception('Error, need at least 2 ids to search for matches in our list of series')

        if use_comparison_df:
            search_df = self.series_comparison_df
        else:
            search_df = self.series_df

        if column_name not in search_df.columns:
            raise Exception('Error, specified column must be present in df')

        # get the set of series_ids that contain the first id
        # self.logger.debug("Searching for: %s" % str(id_list[0]))
        series_ids = search_df[search_df[column_name] == id_list[0]]['SERIES_ID'].tolist()

        if len(series_ids) == 0:
            return pd.DataFrame({})

        selected_series_df = search_df[search_df['SERIES_ID'].isin(series_ids)]

        # iterate over the remaining fragments and keep series_id only if every fragment present
        for idx, item_id in enumerate(id_list):

            if idx > 0:

                # self.logger.debug("Searching for: %s" % str(id_list[idx]))
                series_ids = search_df[search_df[column_name] == id_list[idx]]['SERIES_ID'].tolist()

                if len(series_ids) == 0:
                    return pd.DataFrame({})

                selected_series_df = selected_series_df[selected_series_df['SERIES_ID'].isin(series_ids)]

        # these are specifically not matched to the input fragment series order as
        # some methods use the info about series of alternate ordering to derive observation stats
        if strict_ordering:

            series_ids = selected_series_df.SERIES_ID.unique()
            for series_id in series_ids:

                # get array of fragment ids
                series_df = search_df[search_df['SERIES_ID'] == series_id]
                # print series_df
                # print "My ID list: ", id_list
                # print "Searching column: ", column_name
                items_idx_loc = [series_df[series_df[column_name] == x].index[0] for x in id_list]

                last_item_idx = 0
                for item_idx in items_idx_loc:
                    if item_idx < last_item_idx:
                        # this series is not in requested order so remove from df
                        self.logger.debug("Removing series %d as it is not ordered correctly" % series_id)
                        selected_series_df = selected_series_df[selected_series_df['SERIES_ID'] != series_id]
                        break

                    last_item_idx = item_idx

        # print selected_series_df
        return selected_series_df

    def search_for_mms_to_extend_molids(self, molid_list, strict_order=True, use_comparison_df=False):
        """
        Main method used to search a pre-generated pandas data table of matched series data for any
        series that can be applied as extensions to any/all series found for the molid_list. This method
        mainly chains together a set of other methods
        :param molid_list:
        :return:
        """
        self.logger.info('Searching for series to extend molids: %s' % str(molid_list))

        list_len = len(molid_list)
        if self.return_series_max_len is None:
            self.return_series_max_len = deepcopy(list_len) + 1

        elif list_len >= self.return_series_max_len:
            self.return_series_max_len = deepcopy(list_len) + 1

        if use_comparison_df is False:
            if self.series_df.shape[0] < 4:
                raise Exception("Error, no series found in series_df. Please generate series before trying to "
                                "search them")
        else:
            if self.series_comparison_df.shape[0] < 4:
                raise Exception("Error, no series found in comparison_series_df. Please check input series files "
                                "before trying to search them")

        # TODO: This will only work if the query id's trigger the creation of a series (a matched pair will not work)
        search_series_df = self.return_series_matching_idlist(molid_list, 'MOL_ID',
                                                              strict_ordering=strict_order,
                                                              # False because we start by extracting series for
                                                              # input mols by querying self.series_df not ref set
                                                              use_comparison_df=False)

        if search_series_df.empty or search_series_df.shape[0] < 3:
            raise Exception("Error, no series were found for query molids in target search series_df.")

        # return dataframe
        if use_comparison_df:
            results_df = pd.DataFrame(dtype='int',
                                      columns=['RESULT_ID', 'RESULT_SERIES_SEQ_ID', 'RESULT_MOL_ID',
                                               'RESULT_CONTEXT_ID', 'RESULT_FRAG_ID', 'QUERY_FRAG_ID',
                                               'QUERY_MOL_ID', 'RESULT_MOL_ACTIVITY', 'SOURCE_FILE'])
        else:
            results_df = pd.DataFrame(dtype='int',
                                      columns=['RESULT_ID', 'RESULT_SERIES_SEQ_ID', 'RESULT_MOL_ID',
                                               'RESULT_CONTEXT_ID', 'RESULT_FRAG_ID', 'QUERY_FRAG_ID',
                                               'QUERY_MOL_ID', 'RESULT_MOL_ACTIVITY'])

        result_id = 0

        # now extract the ordered fragment lists for these series and use each one as a query to find matching series
        search_series_lst = search_series_df.SERIES_ID.unique()
        self.logger.info("Found these series ids to use as search seed: %s" % str(search_series_lst))
        for search_series_id in search_series_lst:

            ###############################
            # Convert each series id into a list of fragments we can use to query for similar series
            # 1. get all data for the series as df
            # 2. get fragid for each of the molids
            temp_df = search_series_df[search_series_df.SERIES_ID == search_series_id]

            search_series_context_id = temp_df['CONTEXT_ID'].max()
            self.logger.debug("Query series has context: %s" % search_series_context_id)

            search_frag_ids = []
            for molid in molid_list:
                search_frag_ids.extend(temp_df[temp_df.MOL_ID == molid]['FRAG_ID'].tolist())

            self.logger.info("For input mols found series %s" % str([self.refsmi_dict[x] for x in search_frag_ids]))
            self.logger.info("Will now use this to search for longer series from other mols")

            ###############################
            # Get a list of the series that match our query fragids
            matching_series_df = self.return_series_matching_idlist(search_frag_ids, 'FRAG_ID',
                                                                    strict_ordering=strict_order,
                                                                    use_comparison_df=use_comparison_df)

            # print search_frag_ids, str([self.refsmi_dict[x] for x in search_frag_ids])

            if matching_series_df.empty is True:
                continue

            ################################################
            #  debug to help generate test data
            ################################################
            # result_df_print = copy.deepcopy(matching_series_df)
            # print result_df_print.to_csv()
            # result_df_print['CONTEXT_ID'] = result_df_print['CONTEXT_ID'].apply(lambda x: self.refsmi_dict[x])
            # result_df_print['FRAG_ID'] = result_df_print['FRAG_ID'].apply(lambda x: self.refsmi_dict[x])
            # result_df_print.rename(columns={'CONTEXT_ID': 'CONTEXT',
            #                                 'FRAG_ID': 'FRAG'},
            #                                inplace=True)
            # print search_frag_ids
            # print result_df_print.to_csv()
            #################################################
            # Done debug
            #################################################

            matching_series_lst = matching_series_df.SERIES_ID.unique()
            self.logger.debug("Found %d series that match query series" % len(matching_series_lst))

            ###############################
            # For each series that has matched our query fragid list:
            # 1. Drop it if it's the same series id as our original query
            # 2. Check it has a length longer than our query and last element is new i.e.: extends query
            # 3. Merge results into results dataframe with additional columns for query and context ids

            for matching_series_id in matching_series_lst:

                # make sure our matching series is not the same as one of our query ones
                # only in the case where we are working on the same df that the query came from
                if use_comparison_df is False and matching_series_id in search_series_lst:
                    self.logger.debug("Dropped this series as it is part of the query series set")

                else:
                    self.logger.debug("For molid list %s I have query series fragids %s" %
                                      (str(molid_list), str(search_frag_ids)))

                    matching_series_fragids = matching_series_df[matching_series_df.SERIES_ID ==
                                                                 matching_series_id]['FRAG_ID'].tolist()

                    # check this matched series can be used to extend our existing series
                    if (len(matching_series_fragids) > len(molid_list)) \
                            and (matching_series_fragids[-1] not in search_frag_ids):
                        # this line can be used if the series are always ordered
                        #    and (matching_series_fragids[-1] != search_frag_ids[-1]):
                        # this line in an alternate version for non-ordered series
                        #    and (matching_series_fragids[-1] is not in search_frag_ids):

                        ################################
                        #  these are the results we will pivot into a return df
                        matching_series_molids = matching_series_df[matching_series_df.SERIES_ID ==
                                                                    matching_series_id]['MOL_ID'].tolist()
                        matching_series_seqids = matching_series_df[matching_series_df.SERIES_ID ==
                                                                    matching_series_id]['SERIES_SEQ_ID'].tolist()
                        matching_series_act = matching_series_df[matching_series_df.SERIES_ID ==
                                                                 matching_series_id]['ACTIVITY'].tolist()
                        matching_series_frag_smi = [self.refsmi_dict[x] for x in matching_series_fragids]
                        context_id = matching_series_df[matching_series_df.SERIES_ID ==
                                                        matching_series_id]['CONTEXT_ID'].unique()[0]
                        num_items = len(matching_series_molids)
                        self.logger.debug("The common context SMI is: %s" % self.refsmi_dict[context_id])

                        result_id += 1

                        ################################
                        # these are also added to the pandas table and flag which items in the matching series
                        # were part of the query and which were not (NaN). It's not that intelligent and might
                        # flag an item twice
                        matching_series_query_molids = []
                        matching_series_query_fragids = []

                        for idx, item in enumerate(matching_series_fragids):
                            try:
                                locn = search_frag_ids.index(item)
                                matching_series_query_fragids.append(item)
                                matching_series_query_molids.append(molid_list[locn])
                            except:
                                matching_series_query_fragids.append(0)
                                matching_series_query_molids.append(0)

                        ################################
                        # finally form a return dataframe from all the lists
                        a_result = pd.DataFrame({'RESULT_ID': [result_id] * num_items,
                                                 'RESULT_SERIES_SEQ_ID': matching_series_seqids,
                                                 'RESULT_MOL_ID': matching_series_molids,
                                                 'RESULT_CONTEXT_ID': [context_id] * num_items,
                                                 'QUERY_CONTEXT_ID': [search_series_context_id] * num_items,
                                                 'RESULT_FRAG_ID': matching_series_fragids,
                                                 'QUERY_FRAG_ID': matching_series_query_fragids,
                                                 'QUERY_MOL_ID': matching_series_query_molids,
                                                 'RESULT_MOL_ACTIVITY': matching_series_act
                                                 },
                                                columns=['RESULT_ID', 'RESULT_SERIES_SEQ_ID', 'RESULT_MOL_ID',
                                                         'RESULT_CONTEXT_ID', 'RESULT_FRAG_ID', 'QUERY_FRAG_ID',
                                                         'QUERY_MOL_ID', 'QUERY_CONTEXT_ID', 'RESULT_MOL_ACTIVITY'],
                                                # does not seem to work
                                                # dtype='int'
                                                )

                        if use_comparison_df:
                            a_result['SOURCE_FILE'] = matching_series_df[matching_series_df.SERIES_ID ==
                                                                         matching_series_id]['SOURCE_FILE'].unique()[0]

                        # py2>3 explicit sort=False https://github.com/pandas-dev/pandas/issues/4588#issue-18183895
                        # use False because I control column order and should be aligned
                        results_df = results_df.append(a_result, sort=False)

                        # sanity check
                        # for idx, item in enumerate(matching_series_fragids):
                        #    print result_id, idx+1, matching_series_molids[idx], context_id, item,
                        #    if item in search_frag_ids:
                        #        print item, matching_series_molids[idx]
                        #    else:
                        #        print NaN, NaN
                        # different format sanity check
                        self.logger.debug("Found a series we can extend your query series with:")
                        self.logger.debug("Query (Mols): %s" % str(molid_list))
                        self.logger.debug("Query (Frags): %s" % str([self.refsmi_dict[x] for x in search_frag_ids]))
                        self.logger.debug("Match (Frags): %s" % matching_series_frag_smi)
                        self.logger.debug("Match (Mols): %s" % matching_series_molids)
                        # print ("Found a series we can extend your query series with:")
                        # print ("Query (Mols): %s" % str(molid_list))
                        # print ("Query (Frags): %s" % str([self.refsmi_dict[x] for x in search_frag_ids]))
                        # print ("Match (Frags): %s" % matching_series_frag_smi)
                        # print ("Match (Mols): %s" % matching_series_molids)

                    else:
                        self.logger.debug("Dropped this series as it is the same length as our query")

        if results_df.empty is False:
            results_df = results_df.astype({'RESULT_ID': int,
                                            'RESULT_SERIES_SEQ_ID': int,
                                            'RESULT_MOL_ID': int,
                                            'QUERY_MOL_ID': int,
                                            'RESULT_CONTEXT_ID': int,
                                            'QUERY_CONTEXT_ID': int,
                                            'RESULT_FRAG_ID': int,
                                            'QUERY_FRAG_ID': int
                                            })

            # this is not essential but if converting from df to dict you lose data
            # due to non-unique index, which in turn caused bugs in the unittest
            # results_df['QUERY_MOL_ID'] = results_df[['QUERY_MOL_ID']].astype(int, )
            results_df.reset_index(drop=True, inplace=True)
            # force NaN to zero's. Many later methods use filtering of column values >0 to determine the existence of
            # a result value (e.g.: find all rows with a query frag/mol). Also, NaN is not valid for float columns in pd
            # use zero as surrogate NaN.
            results_df.fillna(value=0, inplace=True)

        return results_df

    ############################################################################################
    #
    # Scoring Methods
    #
    ############################################################################################

    def _iterator_series_pstats(self, frag_list, series_dataframe, strict_order=False):
        """
        Score all examples of frag_list (a series) found in series_dataframe using the p-value
        method of OBoyle, Bostron, Sayle, Gill JMC 2014
        :param frag_list: list of fragments in
        :param series_dataframe: series_dataframe containing all series involving frag_list
        :param strict_order: Typically False as you want to get back the series stats for all
         possible orderings of the frag_list.  If set to true the stats will not change, all
         possible orderings of the frag_list will still be explored but only the specific order
         requested in the frag_list will be returned (i.e.: less return data only)
        :return: a generator containing each series and it's stats
        """
        self.logger.info("searching for: %s and generating stats" % str(frag_list))

        ###########################################################################################
        # For the input series there are many possible alternate series of a different ordering
        # but only some of these are observed.  All alternate permutations can be see here:
        # print [x for x in permutations(frag_list)]
        # but what we really want is observed so we have to transform the series_dataframe into a
        # dict of an_observed_series => {series_id: series_data, series_id: series_data}
        # From this we can determine the series stats as per OBoyle, Bostron, Sayle, Gill JMC 2014
        ###########################################################################################
        observed_series = {}

        # bug: for some reason the .values method returns float values from an int columns
        # which means later id matches are comparing int and float and fail to match!!
        # added [int(y) for y in x] as we get back list of float not int

        # TODO: This is a hack due to a bug in way we score series, remove it and code will fail
        # if series_dataframe.empty:
        #   this series has never been seen!
        #   typically because we took a fragment from source set and transferred onto target set
        #   series, observations, total_series_obs, observed_prob, >expected_prob, enrichment, p_value,
        #   p_value_corrected, significant

        #   self.logger.debug('Series has never been observed!')
        #   yield frag_list, 0, 0, 0, 1.00/float(factorial(len(frag_list))), 0, 0, 0, 0

        mylist = [tuple([int(y) for y in x]) for x in
                  series_dataframe[series_dataframe['FRAG_ID'].isin(frag_list)]
                  [['SERIES_ID', 'FRAG_ID', 'SERIES_SEQ_ID']].values]

        # transpose to dict keyed by series with values eq to list of tuples (frag_id, series_seq_id)
        all_series = {}
        for k, g in groupby(mylist, itemgetter(0)):
            all_series[k] = [el[1:] for el in g]

        # print all_series

        # sort each series by series_seq_id, transpose to dict keyed by sorted (observed!) series
        for series_id in all_series:
            # sort by series_seq_id, return only fragid
            # http://stackoverflow.com/questions/10695139/sort-a-list-of-tuples-by-2nd-item-integer-value
            series_sorted = sorted(all_series[series_id], key=lambda lst: lst[1])
            series_sorted = tuple(map((lambda tup: tup[0]), series_sorted))

            self.logger.debug("Series id %s is sorted to %s" % (series_id, str(series_sorted)))

            # store
            if series_sorted not in observed_series:
                observed_series[series_sorted] = {}

            observed_series[series_sorted][series_id] = None

        ###################################################################
        # Generate the statistics on how significant this series is
        ###################################################################

        total_series_obs = len(all_series)
        n_factorial = factorial(len(frag_list))
        # get p-value for specific observations
        for series in observed_series:

            # enrichment
            observations = len(observed_series[series])
            observed_prob = observations/float(total_series_obs)

            # added try/except due to long series total_series_obs >= 270 causing Long Int to float conversion fail
            # https://stackoverflow.com/questions/13978556/
            # python-long-int-too-large-to-convert-to-float-when-calculating-pi
            #series_too_large = False
            try:
                expected_prob = 1.00/n_factorial
            except Exception:
                yield series, observations, total_series_obs, observed_prob, \
                      np.nan, np.nan, np.nan, np.nan, np.nan
                # yield does not imply continue but we're done so skip to next iteration
                continue

            enrichment = observed_prob / expected_prob

            # p-value (corrected)
            # use binom_test(count, nTrials, pEvent) where pEvent is the expected prob (1/n_factorial)
            p_value = binom_test(observations, total_series_obs, expected_prob)
            # Bonferroni correction, multiplying the original p-value by the degrees of freedom
            # i.e.: p_value * (N! - 1)
            p_value_corrected = p_value * (n_factorial - 1)

            if p_value_corrected <= 0.05:
                significant = True
            else:
                significant = False

            if strict_order:
                if series == tuple(frag_list):
                    yield series, observations, total_series_obs, observed_prob, expected_prob, enrichment, \
                          p_value, p_value_corrected, significant
                else:
                    self.logger.debug("Removed series as non-matching order %s versus %s" % (str(series),
                                                                                             str(tuple(frag_list))))
            else:
                yield series, observations, total_series_obs, observed_prob, expected_prob, enrichment, \
                      p_value, p_value_corrected, significant

    def _iterator_series_scoringmetrics(self, mol_id_list, selected_series_df, apply_filter=False):
        """
        Method will iterate over the df and take each pair of series, convert MolID's to activity values
        then use the ordered series to calculate an cRMSD (centered RMSD). Originally designed to take
        the output from the method search_for_mms_to_extend_molids and process to get cRMSD values.
        Scoring (JAL):
         cRMSD score is from Matched Molecular Series: Measuring SAR Transferability - Christian Kramer,
         Emanuel Ehmki, (Roche), Seventh Joint Sheffield Conference on Chemoinformatics 4th July 2016 and
         J Chem Inf Model. 2017 May 1. doi: 10.1021/acs.jcim.6b00709
         http://cisrg.shef.ac.uk/shef2016/conference-programme/poster-session/...
         matched-molecular-series-measuring-sar-transferability/
         http://pubs.acs.org/doi/abs/10.1021/acs.jcim.6b00709
        Filtering
         apply_filter approach was taken from Keefer and Chang MedChemComm 2017 to remove "candidate series"
         with poor characteristics, pActivity < 5 or skew < 3
        :param: mol_id_list: list of the molecules used in query, order specific, will be sanity checked
        :param selected_series_df: series_dataframe containing all series involving frag_list
        :param apply_filter: boolean, remove candidate series with poor characteristics
        :return: a generator containing each series and it's cRMSD
        """
        # sanity check the data table for the columns we need
        required_col = ['RESULT_ID', 'QUERY_MOL_ID', 'RESULT_MOL_ID', 'QUERY_FRAG_ID', 'RESULT_FRAG_ID']
        for col in required_col:
            if col not in selected_series_df.columns:
                raise Exception("Input data table does not have required columns: %s" % col)

        # sanity check, the series we are scoring must contain all mols in the query series
        # want the RESULT_FRAG_ID where the QUERY_FRAG_ID is not null
        first_series_df = selected_series_df[selected_series_df['RESULT_ID'] == selected_series_df['RESULT_ID'].min()]
        df_mol_ids = first_series_df[first_series_df['QUERY_MOL_ID'] != 0]['QUERY_MOL_ID'].tolist()

        if sorted(df_mol_ids) != sorted(mol_id_list):
            raise Exception("The matching result mol ids in the data frame do not match the mol ids in mol_id_list")

        ##########################
        # Calculation of delta p value
        # (query) - do this here to avoid repetition inside loop
        query_mols_data = []

        # if self.series_from_file == True:
        #    print self.series_comparison_df.columns
        #    print "Using Comparison >>>"
        #    print [self.series_comparison_df[self.series_comparison_df['MOL_ID'] == x]['ACTIVITY'].max() for x in
        #           mol_id_list]
        #    print "--"
        #    print [self.series_comparison_df[self.series_comparison_df['MOL_ID'] == x]['ACTIVITY'].to_csv() for x in
        #           mol_id_list]
        #    print "--"
        #    print [self.series_comparison_df[self.series_comparison_df['MOL_ID'] == x].to_csv() for x in mol_id_list]
        #    print "--"
        #    self.series_comparison_df[self.series_comparison_df['MOL_ID'] == x]['ACTIVITY']
        #    for x_ in mol_id_list:
        #        query_mols_data.append(self.series_comparison_df[self.series_comparison_df['MOL_ID'] ==
        #           int(x_)]['ACTIVITY'].max())
        # else:
        for x_ in mol_id_list:
            query_mols_data.append(self.mol_data_dict[int(x_)][self.data_column_position])

        query_avg = sum([x for x in query_mols_data]) / float(len(query_mols_data))

        sort_order = {}
        for idx, id_ in enumerate(mol_id_list):
            sort_order[id_] = idx

        # iterate over each of the pairs of query:result matched series
        for result_id in selected_series_df['RESULT_ID'].unique():

            just_this_series_df = selected_series_df[selected_series_df['RESULT_ID'] == result_id]

            new_rgroups_df = just_this_series_df[just_this_series_df['QUERY_MOL_ID'] == 0]
            just_this_series_df = just_this_series_df[just_this_series_df['QUERY_MOL_ID'] > 0]

            # sanity check
            #print("Try this one")
            #print(sorted(set(just_this_series_df['QUERY_MOL_ID'].tolist())), sorted(set(mol_id_list)))
            # converted to set() because of situations where we have chiral enantiomers of mol with same fragment
            # creating a series with the same fragment twice, with matching query mol then repeated twice,
            #
            just_query_mol_ids = just_this_series_df['QUERY_MOL_ID'].tolist()
            if sorted(set(just_query_mol_ids)) != sorted(set(mol_id_list)):
                raise Exception("Data error: input dataframe query mols do not match your input mol_id_list")

            # sort according to the input mol_id_list order
            # TODO: bug - ValueError: Categorical categories must be unique, add filter to remove for now
            # case occurs where multiple identical Mol_IDS exist in column 'QUERY_MOL_ID
            if len(set(just_query_mol_ids)) != len(just_query_mol_ids):
                #print("Ditch this set as multiple identical QUERY_MOL_IDs found")
                #print(sorted(set(just_this_series_df['QUERY_MOL_ID'].tolist())), sorted(set(mol_id_list)))
                continue

            just_this_series_df['QUERY_MOL_ID_CAT'] = pd.Categorical(just_this_series_df['QUERY_MOL_ID'],
                                                                     categories=mol_id_list,
                                                                     ordered=True)

            just_this_series_df = just_this_series_df.sort_values('QUERY_MOL_ID_CAT')
            # would this work instead?
            # just_this_series_df['QUERY_MOL_ID_CAT'] = just_this_series_df.reindex(mol_id_list)

            result_mol_ids = just_this_series_df['RESULT_MOL_ID'].tolist()
            result_frag_ids = just_this_series_df['RESULT_FRAG_ID'].tolist()
            result_mols_data = just_this_series_df['RESULT_MOL_ACTIVITY'].tolist()

            if apply_filter:

                act_data_nparr = np.array(result_mols_data)

                # range in activity
                # print "- range = ", (abs(act_data_nparr.max() - act_data_nparr.min()))
                if (abs(act_data_nparr.max() - act_data_nparr.min())) < 0.5:
                    continue

                # skew, via pearsons second coefficient (alternate method avoiding mode)
                # print "- skew = ", (3 * (act_data_nparr.mean() - np.median(act_data_nparr)))/act_data_nparr.std()
                if (3 * (act_data_nparr.mean() - np.median(act_data_nparr)))/act_data_nparr.std() > 3:
                    continue

            ##########################
            # Calculation of cRMSD / RMSD
            #
            # print query_mols_data, result_mols_data
            # convert arrays to np.arrays
            query_mols_data = np.array(query_mols_data)
            result_mols_data = np.array(result_mols_data)

            result_avg = sum([x for x in result_mols_data]) / float(len(result_mols_data))

            # center the data to compensate for the different core, assay or interaction
            # of the compound with the target (see ref)
            query_mols_data_centered = np.array([(x - query_avg) for x in query_mols_data])
            result_mols_data_centered = np.array([(x - result_avg) for x in result_mols_data])
            #
            rmsd = np.sqrt((query_mols_data - result_mols_data) ** 2).mean()
            rmsd_c = np.sqrt((query_mols_data_centered - result_mols_data_centered) ** 2).mean()

            ##########################
            # Predicted activity:
            # Calculation of p values
            # (difference to the mean or diffMean)
            # Calculation via Linear Regression
            #
            new_fragids = new_rgroups_df['RESULT_FRAG_ID'].tolist()
            # diffMean model method
            new_predict_act_dM = []
            # linear model method
            lm_x = linregress(query_mols_data, result_mols_data)
            new_predict_act_lm = []
            new_predict_lm_rsqd = lm_x.rvalue
            # print lm_x.slope, lm_x.intercept

            for molid in new_rgroups_df['RESULT_MOL_ID'].tolist():

                if self.series_from_file:
                    molid_activity = self.series_comparison_df[self.series_comparison_df['MOL_ID'] ==
                                                               molid]['ACTIVITY'].max()

                else:
                    molid_activity = self.mol_data_dict[molid][self.data_column_position]

                new_predict_act_dM.append(query_avg + (molid_activity - result_avg))
                new_predict_act_lm.append((lm_x.slope * molid_activity) + lm_x.intercept)

            ##########################
            # Distance metrics
            # https://docs.scipy.org/doc/scipy/reference/spatial.distance.html
            # from scipy.spatial.distance import cityblock
            #
            manhattan_dist = cityblock(query_mols_data, result_mols_data)
            manhattan_dist_c = cityblock(query_mols_data_centered, result_mols_data_centered)

            ##########################
            # Calculation of r2
            # gets rank order coeff and p-val as list
            # from scipy.stats import pearsonr, spearmanr
            #
            pearson_r = pearsonr(query_mols_data, result_mols_data)
            spearman_r = spearmanr(query_mols_data, result_mols_data)

            yield result_id, result_mol_ids, result_frag_ids, new_fragids, \
                  new_predict_act_dM, new_predict_act_lm, new_predict_lm_rsqd, rmsd, rmsd_c, \
                  manhattan_dist, manhattan_dist_c, pearson_r[0], pearson_r[1], spearman_r[0], spearman_r[1]

    ############################################################################################
    #
    # Write results methods
    #
    ############################################################################################

    def write_raw_series_to_file(self, csv_out, apply_pre_filter=False):
        """ This method writes out raw matched series to file
        """

        self.logger.info('Writing raw series to file: %s' % csv_out)

        with open(csv_out, "w") as out_file:

            out_file.write('SERIES_ID,SERIES_SEQ_ID,MOL_ID,CONTEXT,FRAG,ACTIVITY\n')

            for series_id, series_seq_id, ctx_id, \
                frag_id_L, molid_L, mol_act in self._iterator_mmp_series_numeric(apply_pre_filter=apply_pre_filter):

                print_string = str(series_id) + "," + str(series_seq_id) + "," + str(molid_L) + ","
                print_string += self.refsmi_dict[ctx_id] + ","
                print_string += self.refsmi_dict[frag_id_L] + ","
                print_string += str(mol_act) + "\n"

                out_file.write(print_string)

        self.logger.info('Done writing file: %s' % csv_out)

    def return_scored_series_dataframe(self,
                                       mol_id_list,
                                       results_dataframe,
                                       return_dataframe,
                                       append=False,
                                       apply_filter=False):
        """
        Take a dataframe as input and validate as it expects the output from search_for_mms_to_extend_molids. Now
        iterate across the result series and score then return....
        :param mol_id_list the mol id list used to generate the results_dataframe
        :param results_dataframe name of input dataframe to reformat/score
        :param return_dataframe name of dataframe to return, this var must exist but can be set to None
        :param append should I append the reformatted results_dataframe to return_dataframe or is it a new df
        :param apply_filter (Filter)
        :return reformatted, scored dataframe
        """
        # self.logger.info('Writing scored series to file %s' % csv_out)

        # validation of df structure
        required_col = ['RESULT_ID', 'RESULT_SERIES_SEQ_ID', 'QUERY_MOL_ID', 'RESULT_MOL_ID',
                        'QUERY_FRAG_ID', 'QUERY_CONTEXT_ID', 'RESULT_CONTEXT_ID', 'RESULT_FRAG_ID']
        for col in required_col:
            if col not in results_dataframe.columns:
                raise Exception("Input data table does not have required columns: %s" % col)

        inc_source_sol = False
        if 'SOURCE_FILE' in results_dataframe.columns:
            inc_source_sol = True

        # catch for empty table
        if results_dataframe.shape[0] == 0:
            print("No results found")
            return False

        # get RMSE scoring (work on unexpanded series)
        results_container = {}
        unique_series = set()

        for result_id, result_mol_ids, result_frag_ids, new_fragids, \
            new_predict_act_dM, new_predict_act_lm, new_predict_lm_rsqd,\
            rmsd, rmsd_c, manhattan_dist, manhattan_dist_c, \
            pearson_r2, pearson_p, spearman_r2, spearman_p \
                in self._iterator_series_scoringmetrics(mol_id_list, results_dataframe, apply_filter=apply_filter):

            results_container[result_id] = [result_id, result_mol_ids, result_frag_ids, new_fragids,
                                            new_predict_act_dM, new_predict_act_lm, new_predict_lm_rsqd,
                                            rmsd, rmsd_c, manhattan_dist, manhattan_dist_c,
                                            pearson_r2, pearson_p, spearman_r2, spearman_p]
            unique_series.add(result_id)

        # debugging - need to see what I have
        # with open('rmse_score.csv', 'wb') as outrmse:
        #   for key, value in results_container.items():
        #        outrmse.write(str(key)+","+str(value)+"\n")
        #        print key, value

        # get all new expanded series - basically the ideas of things we should make
        # ...and get the p-value based scoring for them
        # removed this as if we do filtering on the result series by range and skew, we might lose result_id's
        # instead build this in the above loop
        # unique_series = results_dataframe['RESULT_ID'].unique()

        #
        # setup return dataframe inc headers
        #
        if append is False:
            return_dataframe = None

        # num_items = len(mol_id_list)
        header_row = ['LAST_QUERY_SERIES_MOLID', 'LAST_QUERY_SERIES_MOLSMI', 'QUERY_MOL_CONTEXT',
                      'RESULT_MOL_CONTEXT', 'QUERY_MOL_FRAG_L', 'NEW_FRAG_R']
        # series_frag_smi list
        # (a) here we only have the query frags
        for idx in range(self.return_series_max_len):
            header_row.append('QUERY_SERIES_FRAG_' + str(idx + 1))
        header_row.extend(['SERIES_N_OBS', 'SERIES_N_TOTALOBS', 'OBS_PROB', 'EXPECTED_OBS', 'ENRICHMENT',
                           'P_VALUE', 'P_VALUE_CORRECTED', 'SIGNIFCANCE'])
        # results_container[result_id]
        # (b) here we have the query frag + 1 additional matched frag that extends the series hence +1
        for idx in range(self.return_series_max_len + 1):
            header_row.append('RESULT_SERIES_MOLID_' + str(idx + 1))
        for idx in range(self.return_series_max_len + 1):
            header_row.append('RESULT_SERIES_FRAGSMI_' + str(idx + 1))
        header_row.extend(['PRED_ACT_DDP', 'PRED_ACT_LR', 'PRED_ACT_LR_R2', 'RMSE', 'cRMSE', 'MD', 'cMD',
                           'PEARSON_R2', 'PEARSON_PVAL', 'SPEARMAN_R2', 'SPEARMAN_PVAL'])
        if inc_source_sol:
            header_row.append('SOURCE_FILE')

        #
        results_rows = []

        # iterate over the results, foreach series get the last molecule from query and the new fragment
        # to make (can be more than one per series), return single new fragment idea per line with score
        for result_id in unique_series:

            result_subdf = results_dataframe[results_dataframe['RESULT_ID'] == result_id]
            query_context_id = result_subdf['QUERY_CONTEXT_ID'].max()
            result_context_id = result_subdf['RESULT_CONTEXT_ID'].max()
            if inc_source_sol:
                source_file = result_subdf.SOURCE_FILE.unique()[0]
                # print source_file

            # get rowid for last element in query
            last_query_frag_seqid = result_subdf[result_subdf['QUERY_FRAG_ID'] != 0]['RESULT_SERIES_SEQ_ID'].max()
            query_fragids = result_subdf[result_subdf['QUERY_FRAG_ID'] != 0]['QUERY_FRAG_ID'].tolist()
            query_fragids = [int(x) for x in query_fragids]
            last_query_fragid = query_fragids[-1]
            last_query_molid = result_subdf[result_subdf['QUERY_FRAG_ID'] != 0]['QUERY_MOL_ID'].tolist()[-1]
            # get all result frag_id, mol_id that come after the last element in query
            new_fragments = [tuple(x) for x in result_subdf[result_subdf['RESULT_SERIES_SEQ_ID'] >
                                                            last_query_frag_seqid][['RESULT_FRAG_ID',
                                                                                    'RESULT_MOL_ID']].values]

            # I need to do p-scoring for the frag list, which needs all example series matching frag list.
            # These already exist in the results_dataframe but other series might also exist that were found
            # from a different set of query frags. Therefore we need to filter the results_dataframe to get only
            # all series matching the frag list.
            for (new_fragid, new_molid) in new_fragments:
                result_fragids = deepcopy(query_fragids)
                result_fragids.append(int(new_fragid))

                self.logger.debug('New frag %s searching comparison dict %s' % (new_fragid, self.series_from_file))

                # For the new extended series that we have found (query plus new frag_id), score each one
                # TODO: Look at why I have to use a dataframe level property, inconsistent with use_comparison_df param
                # in other words I should be able to just pass use_comparison_df, maybe make use_comparison_df an object
                # level property
                if self.series_from_file is False:
                    temp_results = self.return_series_matching_idlist(result_fragids, 'FRAG_ID')
                else:
                    temp_results = self.return_series_matching_idlist(result_fragids, 'FRAG_ID', use_comparison_df=True)
                    # TODO: you should not be looking in comparison dict as data came from source
                    if self.running_fullscan and temp_results.empty:
                        temp_results = self.return_series_matching_idlist(result_fragids, 'FRAG_ID')

                # TODO: if using non strict matching, this makes no sense, return nulls instead
                for series, n_obs, tot_series_obs, obs_prob, exptd_prob, enrichment, p_val, \
                    p_val_corrected, signif in self._iterator_series_pstats(result_fragids, temp_results,
                                                                            strict_order=True):
                    #
                    single_row = []
                    frag_l = self.refsmi_dict[last_query_fragid]
                    frag_r = self.refsmi_dict[new_fragid]
                    query_context_smi = self.refsmi_dict[query_context_id]
                    result_context_smi = self.refsmi_dict[result_context_id]
                    series_frag_smi = [self.refsmi_dict[x] for x in series]
                    series_frag_smi = series_frag_smi[:-1]
                    this_series_len = len(series_frag_smi)
                    series_frag_smi.extend([None for _ in range(this_series_len, self.return_series_max_len)])

                    # TODO: This should be the smiles not the context
                    if self.series_from_file:
                        last_query_mol_smi = self.series_comparison_df[self.series_comparison_df['MOL_ID'] ==
                                                                       last_query_molid]['CONTEXT_ID'].max()
                    else:
                        last_query_mol_smi = self.mol_smi_dict[last_query_molid]

                    # print last_query_molid, last_query_mol_smi, frag_l, "-->", frag_r, series_frag_smi, n_obs, \
                    #    tot_series_obs, obs_prob, exptd_prob, enrichment, p_val, \
                    #    p_val_corrected, signif, results_container[result_id]
                    single_row.extend([last_query_molid, last_query_mol_smi, query_context_smi,
                                       result_context_smi, frag_l, frag_r])
                    single_row.extend(series_frag_smi)
                    #
                    single_row.extend([n_obs, tot_series_obs, obs_prob, exptd_prob, enrichment,
                                       p_val, p_val_corrected, signif])
                    # get the result molids that match the query molids via fragment series match
                    # and add the new mol that adds new information
                    single_row.extend(results_container[result_id][1])
                    single_row.append(new_molid)
                    single_row.extend([None for _ in range(this_series_len, self.return_series_max_len)])
                    single_row.extend([self.refsmi_dict[x] for x in results_container[result_id][2]])
                    single_row.append(frag_r)
                    single_row.extend([None for _ in range(this_series_len, self.return_series_max_len)])
                    # get the predicted activity for this fragment
                    idx_loc = results_container[result_id][3].index(new_fragid)
                    single_row.append(results_container[result_id][4][idx_loc])
                    single_row.append(results_container[result_id][5][idx_loc])
                    # same again for the fragments of this series
                    self.logger.debug("frag id %s has activity %s" % (new_fragid,
                                                                      results_container[result_id][4][idx_loc]))
                    # rest of the info
                    single_row.extend(results_container[result_id][6:])
                    if inc_source_sol:
                        single_row.append(source_file)
                    results_rows.append(single_row)
                    # print result_id, last_query_molid, last_query_fragid, "-->", new_fragid, series, n_obs, \
                    #    tot_series_obs, obs_prob, exptd_prob, enrichment, p_val, \
                    #    p_val_corrected, signif, results_container[result_id]

        if append is False:
            return_dataframe = pd.DataFrame(results_rows, columns=header_row)

        else:
            return_dataframe = return_dataframe.append(pd.DataFrame(results_rows, columns=header_row), sort=False)
            # return_dataframe = return_dataframe.rename_axis(None)

        # add enumerated product
        #print(return_dataframe.shape)
        if return_dataframe.shape[0] < 1:
            return pd.DataFrame(columns=header_row)
        else:
            return_dataframe = self.enumerate_products(return_dataframe, 'QUERY_MOL_CONTEXT', 'NEW_FRAG_R')

        # add max activity for given context
        def get_series_max_act(row):
            return self.max_series_potency[self.refsmi_dict[row['QUERY_MOL_CONTEXT']]]

        return_dataframe['MAX_SERIES_ACT'] = return_dataframe.apply(lambda row: get_series_max_act(row), axis=1)

        # reorder headers
        header_row.extend(['ENUMERATED_PRODUCT', 'NOVEL', 'MAX_SERIES_ACT'])
        # IL replaced deprecated .reindex_axis with .reindex
        return_dataframe = return_dataframe.reindex(header_row, axis=1)
        #
        return return_dataframe

    def write_series_dataframe(self, dataframe, csv_out):
        """
        :param dataframe: input dataframe to write to file
        :param csv_out: name of file to write df to
        Only real use for this right now is for object level access, debugging,
        avoiding pandas import at script wrapper level
        """
        dataframe.to_csv(csv_out, index=False, float_format='%.3f')

    ############################################################################################
    #
    # Methods to parse a set of series files from disk
    #
    ############################################################################################

    def set_base_dir(self, base_dir):
        """Set the base directory for this script to run from, all files stored here"""

        if os.path.isdir(base_dir):

            self.series_base_dir = base_dir

            if self.series_base_dir[-1] != "/":
                self.series_base_dir += "/"

            self.logger.info('Base Directory set as %s' % self.series_base_dir)
            # print ('Base Directory set as %s' % self.series_base_dir)

        else:
            self.logger.warn('Fatal: Base Directory does not exist %s' % self.series_base_dir)
            raise Exception('Error Base Directory does not exist')

    def get_series_filelist_from_dir(self):
        """get list of *.series files from base_dir,
        if present store in crc_raw_file_list & return True"""
        #
        if self.series_base_dir is None:
            self.logger.warn('Fatal: Base Directory not set %s')
            raise Exception('Error Base Directory not set')

        self.series_file_list = glob.glob(self.series_base_dir + '*.series')
        self.logger.info('Found %d files' % len(self.series_file_list))

        if len(self.series_file_list) >= 1:
            return True

        else:
            return False

    def parse_directory_of_series_files(self):
        """
        Read all files in directory with .series extension and merges into one huge table
        :param: dir_name:
        :return: None
        """
        if self.series_base_dir is None or len(self.series_file_list) < 1:
            self.logger.warn('Fatal: Base Directory not set %s')
            raise Exception('Error Base Directory not set')

        self.logger.info('Parsing dir of files from %s' % self.series_base_dir)

        self.ref_series_df = pd.DataFrame([], columns=['SERIES_ID', 'SERIES_SEQ_ID', 'CONTEXT',
                                                       'FRAG', 'MOL_ID', 'ACTIVITY'])

        required_col = ['SERIES_ID', 'SERIES_SEQ_ID', 'CONTEXT', 'FRAG', 'MOL_ID', 'ACTIVITY']
        max_series_id = 0

        for series_file in self.series_file_list:

            # print series_file
            temp_df = pd.read_csv(series_file)  # , index_col=False)
            # print temp_df.columns

            # sanity check the data table for the columns we need
            for col in required_col:
                if col not in temp_df.columns:
                    raise Exception("Input CSV %s does not have required columns: %s" % (series_file, col))

            # re-sequence the series ID's
            if max_series_id == 0:
                max_series_id = temp_df['SERIES_ID'].max()
            else:
                max_series_id = self.ref_series_df['SERIES_ID'].max()
            # print max_series_id

            temp_df['SERIES_ID'] = temp_df['SERIES_ID'] + max_series_id
            temp_df['SOURCE_FILE'] = os.path.basename(series_file)

            # py2>3 explicit sort=False added
            self.ref_series_df = self.ref_series_df.append(temp_df, sort=False)
            self.logger.info('Appended dataframe shape %s to master dataframe %s' %
                             (str(temp_df.shape), str(self.ref_series_df.shape)))
            # print ('Appended dataframe shape %s to master dataframe %s' % (str(temp_df.shape),
            #   str(self.ref_series_df.shape)))
            # print self.ref_series_df['SERIES_ID'].max()

        self.series_comparison_df = self.ref_series_df

    def setup_pregenerated_series_data_for_mms(self, dir_name):
        """
        Runs a set of methods to read in a dir of series files ready for searching
        :param dir_name: a directory of series files, full path
        :return:
        """
        self.logger.info('Setting up searching using pre-generated series data')

        self.set_base_dir(dir_name)
        self.get_series_filelist_from_dir()
        self.parse_directory_of_series_files()
        self.series_from_file = True

        self.logger.info('Parsed data to series data structure of shape %s' % str(self.series_comparison_df.shape))

        def convert_smi_to_id(smi_str):
            """convert a fragment or context smiles to a reference id consistent with existing mmp object ids"""

            if smi_str in self.refsmi_dict:
                frag_id = self.refsmi_dict[smi_str]

            else:
                self.refsmi_id += 1
                self.refsmi_dict[smi_str] = self.refsmi_id
                self.refsmi_dict[self.refsmi_id] = smi_str
                frag_id = self.refsmi_id

            return frag_id

        # flip 'CONTEXT' and 'FRAG' columns to reference ID values
        # print self.series_comparison_df.shape
        # print self.series_comparison_df.head()
        self.series_comparison_df['CONTEXT'] = self.series_comparison_df['CONTEXT'].apply(
            lambda x: convert_smi_to_id(x))
        self.series_comparison_df['FRAG'] = self.series_comparison_df['FRAG'].apply(lambda x: convert_smi_to_id(x))

        # rename them as ID columns as no longer raw smi strings
        self.series_comparison_df.rename(columns={'CONTEXT': 'CONTEXT_ID',
                                                  'FRAG': 'FRAG_ID'},
                                         inplace=True)

        # print self.series_comparison_df.to_csv()

    def enumerate_products(self, dataframe, context_col, fragment_col):
        """
        Uses MMPEnumerateNewMols object to enumerate new molecules from ctx + frag pair
        which are expected to already be labelled with isotopic labelling to highlight
        the connection points
        :param dataframe:
        :param context_col:
        :param fragment_col:
        :return:
        """
        # create dict
        self.enumerated_products_smi = {}
        for (ctx, frag) in zip(dataframe[context_col].tolist(), dataframe[fragment_col].tolist()):
            self.enumerated_products_smi[(ctx, frag)] = None

        # function to return value from this dict:
        def get_product(row):
            return self.enumerated_products_smi[(row[context_col], row[fragment_col])]

        # function to return value for novelty
        def get_novelty(row):
            if row['ENUMERATED_PRODUCT'] in self.mol_smi_dict:
                return self.mol_smi_dict[row['ENUMERATED_PRODUCT']]
            else:
                return True

        # write rxn files
        enum_object = enum_mols.MMPEnumerateNewMols(self.logger)
        enum_object.write_rxn_files()
        enum_object.write_reactants_simple_dict(self.enumerated_products_smi)
        enum_object.do_reactions()

        for cut_type, rtn_ctx, rtn_frag, new_mol in enum_object.yield_products_simple_dict_input():
            if (rtn_ctx, rtn_frag) in self.enumerated_products_smi:
                self.enumerated_products_smi[(rtn_ctx, rtn_frag)] = new_mol
            else:
                self.logger.debug("got a return I was not expecting: %s, %s -> %s" % (rtn_ctx, rtn_frag, new_mol))

        # standardise
        temp_smifi = tempfile.NamedTemporaryFile(suffix=".smi", delete=False, encoding='utf-8', mode='wt')
        std_smi_lookup = {}
        # take the product smi and store as product_smi => None
        arbitary_id = 0
        for smi_no_std in list(self.enumerated_products_smi.values()):
            arbitary_id += 1
            std_smi_lookup[smi_no_std] = arbitary_id
            std_smi_lookup[arbitary_id] = smi_no_std
            temp_smifi.write(smi_no_std + " " + str(arbitary_id) + "\n")
        temp_smifi.close()

        #
        # send in a dict of old_smi => old_smi
        # should turn this into old_smi => new_smi
        self.logger.debug("Requested standardisation on %s" % temp_smifi.name)
        new_smi_dict = self.generate_std_smiles(temp_smifi.name, smi_id_map='id_smi')
        std_smi_lookup.update(new_smi_dict)

        #
        for key, value in list(self.enumerated_products_smi.items()):
            # print value, " >> ", std_smi_lookup[std_smi_lookup[value]]
            self.enumerated_products_smi[key] = std_smi_lookup[std_smi_lookup[value]]

        # add the return enumerated mols to the df
        dataframe['ENUMERATED_PRODUCT'] = dataframe.apply(lambda row: get_product(row), axis=1)
        dataframe['NOVEL'] = dataframe.apply(lambda row: get_novelty(row), axis=1)
        #
        del enum_object
        #
        return dataframe

    def auto_search(self, m, n, strict_ordering=False, use_comparison_df=False):
        """
        Method to automatically search a series dataframe for series that can be extended
        to create molecules of activity the same as or greater than the best molecule in
        the set.  The method will remove m items from every series in the set of molecules. It
        will then decompose this set of m items into every possible subset down to length m but
        maintaining the sequence order.
        :param m: the number of items removed from a given series to create the longest seed
        :param n: the smallest number of items to be used in a query series, decomposed from m
        :return: new idea compounds as iterator
        """
        self.logger.info('Initialised auto_search with parameters m:%s, n:%s' % (m, n))

        if n > m:
            raise Exception("Error, m must be greater or equal to n")
        else:
            self.return_series_max_len = m

        series_ids = self.series_df['SERIES_ID'].unique()
        all_results = pd.DataFrame(dtype='int',
                                   columns=['RESULT_ID', 'SERIES_ID', 'RESULT_SERIES_SEQ_ID', 'RESULT_MOL_ID',
                                            'RESULT_CONTEXT_ID', 'RESULT_FRAG_ID', 'QUERY_CONTEXT_ID', 'QUERY_FRAG_ID',
                                            'QUERY_MOL_ID', 'QUERY_ORDER', 'RESULT_MOL_ACTIVITY'])
        iteration = 0
        result_id = 0

        # iterate over every matched series, use last m values as search query
        for series_id in series_ids:

            series_fragids = self.series_df[self.series_df['SERIES_ID'] == series_id]['FRAG_ID'].tolist()
            series_molids = self.series_df[self.series_df['SERIES_ID'] == series_id]['MOL_ID'].tolist()
            query_context_id = self.series_df[self.series_df['SERIES_ID'] == series_id]['CONTEXT_ID'].max()
            # self.logger.debug("Query context: %s" %
            #                  self.refsmi_dict[self.series_df[self.series_df['SERIES_ID'] == series_id]
            #                  ['CONTEXT_ID'].max()])

            if len(series_fragids) > m:
                series_fragids = series_fragids[-m:]
                series_molids = series_molids[-m:]

            series_frag_n_mols = [tuple([val, series_molids[idx]]) for idx, val in enumerate(series_fragids)]

            # get every length of series we will generate between m and n then generate
            # the actual ordered sub-series from a series using combinations, start with
            # the longest series first as later we will drop/never process a series if it's
            # series_id has been seen before, so take longest series first
            for size in range(m, n - 1, -1):

                for query_tuples in combinations(series_frag_n_mols, size):

                    iteration += 1
                    query_fragids = [x[0] for x in query_tuples]
                    query_molids = [x[1] for x in query_tuples]

                    # run search, remove the current series_id from result list
                    return_df = self.return_series_matching_idlist(query_fragids, 'FRAG_ID',
                                                                   strict_ordering=strict_ordering,
                                                                   use_comparison_df=use_comparison_df)
                    #
                    # search might return nothing, or might only return the series the query came from
                    if not return_df.empty:
                        return_df = return_df[return_df['SERIES_ID'] != series_id]
                    else:
                        # this should never happen for internal SAR transfer but can for external set
                        continue
                    #
                    if not return_df.empty:

                        for sgl_rtn_series_id in return_df['SERIES_ID'].unique():

                            sgl_rtn_return_df = return_df[return_df['SERIES_ID'] == sgl_rtn_series_id]

                            # only add if we've not seen it before
                            # We're processing series in reverse length order so we'll capture the longest series
                            # first and drop anything shorter with the same series_id.
                            #
                            if sgl_rtn_return_df['SERIES_ID'].max() not in all_results['SERIES_ID'].unique():

                                matching_series_fragids = sgl_rtn_return_df['FRAG_ID'].tolist()

                                # check this matched series can be used to extend our existing series
                                if (len(matching_series_fragids) > len(query_fragids)) \
                                    and (matching_series_fragids[-1] not in query_fragids):

                                    result_id += 1
                                    # debugging
                                    # print size, len(query_fragids), len(matching_series_fragids)

                                    # reformat results
                                    matching_series_molids = sgl_rtn_return_df['MOL_ID'].tolist()
                                    matching_series_seqids = sgl_rtn_return_df['SERIES_SEQ_ID'].tolist()
                                    matching_series_act = sgl_rtn_return_df['ACTIVITY'].tolist()
                                    # matching_series_frag_smi = [self.refsmi_dict[x] for x in matching_series_fragids]
                                    matching_series_context_id = sgl_rtn_return_df['CONTEXT_ID'].unique()[0]
                                    matching_series_id = sgl_rtn_return_df['SERIES_ID'].max()
                                    if use_comparison_df:
                                        matching_series_source = sgl_rtn_return_df['SOURCE_FILE'].unique()[0]
                                    num_items = len(matching_series_molids)

                                    # don't know if single or double so skip this
                                    # (ms_ctx1_id, ms_ctx1_id) = inv_cantor(matching_series_context_id)
                                    # self.logger.debug("Result context: %s" %
                                    #                  self.refsmi_dict[matching_series_context_id])

                                    ################################
                                    # these are also added to the pandas table and flag which items in the matching
                                    # series were part of the query and which were not (NaN). It's not that intelligent
                                    # and might flag an item twice
                                    matching_series_query_molids = []
                                    matching_series_query_fragids = []
                                    for idx, item in enumerate(matching_series_fragids):
                                        try:
                                            locn = query_fragids.index(item)
                                            matching_series_query_fragids.append(item)
                                            matching_series_query_molids.append(query_molids[locn])
                                        except:
                                            matching_series_query_fragids.append(0)
                                            matching_series_query_molids.append(0)

                                    # Need the query order added to the dataframe as we did non-ordered matching
                                    matching_series_query_molids_order = []
                                    for idx, query_id in enumerate(matching_series_query_molids):
                                        if query_id in query_molids:
                                            matching_series_query_molids_order.append(query_molids.index(query_id) + 1)
                                        else:
                                            matching_series_query_molids_order.append(0)

                                    # print [result_id] * num_items
                                    # print [matching_series_id] * num_items
                                    # print matching_series_seqids
                                    # print matching_series_molids
                                    # print [matching_series_context_id] * num_items
                                    # print 'matching_series_fragids ', matching_series_fragids
                                    # print 'matching_series_query_fragids', matching_series_query_fragids
                                    # print 'matching_series_query_molids', matching_series_query_molids
                                    # print matching_series_act

                                    ################################
                                    # finally form a return dataframe from all the lists
                                    reformatted_result_df = pd.DataFrame({
                                                             'RESULT_ID': [result_id] * num_items,
                                                             'SERIES_ID': [matching_series_id] * num_items,
                                                             'RESULT_SERIES_SEQ_ID': matching_series_seqids,
                                                             'RESULT_MOL_ID': matching_series_molids,
                                                             'RESULT_CONTEXT_ID': [matching_series_context_id] * num_items,
                                                             'RESULT_FRAG_ID': matching_series_fragids,
                                                             'QUERY_CONTEXT_ID': [query_context_id] * num_items,
                                                             'QUERY_FRAG_ID': matching_series_query_fragids,
                                                             'QUERY_MOL_ID': matching_series_query_molids,
                                                             'QUERY_ORDER': matching_series_query_molids_order,
                                                             'RESULT_MOL_ACTIVITY': matching_series_act
                                                             },
                                                             columns=['RESULT_ID', 'SERIES_ID', 'RESULT_SERIES_SEQ_ID',
                                                                      'RESULT_MOL_ID', 'RESULT_CONTEXT_ID',
                                                                      'RESULT_FRAG_ID', 'QUERY_CONTEXT_ID',
                                                                      'QUERY_FRAG_ID', 'QUERY_MOL_ID',
                                                                      'QUERY_ORDER', 'RESULT_MOL_ACTIVITY']
                                                            )
                                    if use_comparison_df:
                                        reformatted_result_df['SOURCE_FILE'] = matching_series_source
                                        # [matching_series_source] * num_items

                                    all_results = all_results.append(reformatted_result_df)

                                    self.logger.debug("series %s, search iteration %s using frag_id list %s" %
                                                        (series_fragids, series_id, query_fragids))
                                    self.logger.debug("merging result_df (%s) from query %s to all_results (now %s)" %
                                                        (str(return_df.shape), series_fragids, str(all_results.shape)))

                            else:
                                self.logger.debug("Skipped as already have this series in results set")

                    else:
                        self.logger.debug("No results")

        # all_results = all_results.drop_duplicates().reset_index(drop=True)
        # print all_results.columns
        all_results = all_results.astype({'RESULT_ID': int,
                                          'SERIES_ID': int,
                                          'RESULT_SERIES_SEQ_ID': int,
                                          'RESULT_MOL_ID': int,
                                          'RESULT_CONTEXT_ID': int,
                                          'RESULT_FRAG_ID': int,
                                          'QUERY_CONTEXT_ID': int,
                                          'QUERY_FRAG_ID': int,
                                          'QUERY_MOL_ID': int,
                                          'QUERY_ORDER': int
                                         })

        self.logger.info('Auto search done')
        return all_results

    def auto_search_fullscan(self, strict_ordering=False, use_comparison_df=False):
        # TODO: describe this versus auto_search
        """ """

        self.running_fullscan = True
        self.logger.info('Initialised auto_search_fullscan')

        result_id = 0
        all_results = pd.DataFrame(dtype='int',
                                   columns=['RESULT_ID', 'SERIES_ID', 'RESULT_SERIES_SEQ_ID', 'RESULT_MOL_ID',
                                            'RESULT_CONTEXT_ID', 'RESULT_FRAG_ID', 'QUERY_CONTEXT_ID', 'QUERY_FRAG_ID',
                                            'QUERY_MOL_ID', 'QUERY_ORDER', 'RESULT_MOL_ACTIVITY'])

        ############################
        #
        #
        self.logger.info("Getting frag lists for each series")
        series_ids = self.series_df['SERIES_ID'].unique()
        self.logger.info("Got %s series to work on" % len(series_ids))
        series_frags = {}
        target_frags = {}
        for series_id in series_ids:
            # TODO: storing list as it is ordered, but set here would be faster to avoid conversion later?
            series_frags[series_id] = self.series_df[self.series_df['SERIES_ID'] == series_id]['FRAG_ID'].tolist()

        ############################
        #
        #

        # useful counters and things needed in loop:
        total_num_series = len(series_ids)
        trunc = total_num_series
        self.return_series_max_len = self.min_series_length

        self.logger.info("Working on %s series, %s comparisons" % (total_num_series,
                                                                   (total_num_series * total_num_series / 2)))

        if use_comparison_df:
            target_series_ids = self.series_comparison_df['SERIES_ID'].unique()
            target_series_df = self.series_comparison_df
            for series_id in target_series_ids:
                # TODO: storing list as it is ordered, but set here would be faster to avoid conversion later?
                target_frags[series_id] = target_series_df[target_series_df['SERIES_ID'] ==
                                                           series_id]['FRAG_ID'].tolist()
        else:
            target_series_ids = series_ids[-1 * trunc:]
            target_series_df = self.series_df
            target_frags = series_frags

        for series_id_a in series_ids:
            trunc -= 1
            intersect_frags = []

            # TODO: last comparison is last item versus self
            for series_id_b in target_series_ids:

                # series_frags[series_id_a]
                # series_frags[series_id_b]
                interset_frags = set(series_frags[series_id_a]).intersection(set(target_frags[series_id_b]))
                num_interset_frags = len(interset_frags)

                #
                # This is an MMS pair, so now we need to process it
                #
                if num_interset_frags > self.min_series_length:

                    if num_interset_frags > self.return_series_max_len:
                        self.return_series_max_len = deepcopy(num_interset_frags)

                    # TODO: still have cases here where two series are identical and intersection same, so no extension
                    # print ">>>"
                    # print "series: ", series_id_a, " has items ", series_frags[series_id_a]
                    # print "series: ", series_id_b, " has items ", series_frags[series_id_b]
                    # print "with intersect size:", len(interset_frags), "and items ", interset_frags

                    # ordered list of fragments for the series
                    series_a_df = self.series_df[self.series_df['SERIES_ID'] == series_id_a]
                    series_b_df = target_series_df[target_series_df['SERIES_ID'] == series_id_b]

                    series_a_fragids = series_a_df['FRAG_ID'].tolist()
                    series_b_fragids = series_b_df['FRAG_ID'].tolist()
                    series_a_molids = series_a_df['MOL_ID'].tolist()
                    series_b_molids = series_b_df['MOL_ID'].tolist()

                    # the intersection fragments, but ordered by the original series order
                    series_a_query = [x for x in series_a_fragids if x in interset_frags]
                    series_b_query = [x for x in series_b_fragids if x in interset_frags]

                    series_a_df.loc[:, 'QUERY_FRAG_ID'] = [x if x in interset_frags else 0 for x in series_a_fragids]
                    series_b_df.loc[:, 'QUERY_FRAG_ID'] = [x if x in interset_frags else 0 for x in series_b_fragids]

                    # indexes are zero indexed and we use them as a pseudo for sequence order so increment by 1
                    series_a_df.loc[:, 'QUERY_SERIES_SEQID'] = [series_b_query.index(x) + 1 if x in interset_frags else
                                                                0 for x in series_a_fragids]
                    series_b_df.loc[:, 'QUERY_SERIES_SEQID'] = [series_a_query.index(x) + 1 if x in interset_frags else
                                                                0 for x in series_b_fragids]

                    series_a_df['QUERY_MOL_ID'] = [series_b_molids[series_b_query.index(x)] if x in interset_frags else
                                                   0 for x in series_a_fragids]
                    series_b_df['QUERY_MOL_ID'] = [series_a_molids[series_a_query.index(x)] if x in interset_frags else
                                                   0 for x in series_b_fragids]

                    # ensure we can extend the query series
                    # if so then reformat the df
                    if series_a_df['QUERY_SERIES_SEQID'].tolist()[-1] == 0:

                        if use_comparison_df:
                            continue
                        else:
                            result_id += 1
                            series_a_df['RESULT_ID'] = [result_id] * len(series_a_fragids)
                            series_a_df['QUERY_CONTEXT_ID'] = [series_b_df.loc[:, 'CONTEXT_ID'].max()] * \
                                len(series_a_fragids)
                            # named new_df because an in place rename will affect second below if statement,
                            # would no longer find CONTEXT_ID column in series_a_df
                            new_df = series_a_df.rename(columns={'SERIES_SEQ_ID': 'RESULT_SERIES_SEQ_ID',
                                                                 'CONTEXT_ID': 'RESULT_CONTEXT_ID',
                                                                 'FRAG_ID': 'RESULT_FRAG_ID',
                                                                 'MOL_ID': 'RESULT_MOL_ID',
                                                                 'ACTIVITY': 'RESULT_MOL_ACTIVITY',
                                                                 'QUERY_SERIES_SEQID': 'QUERY_ORDER'})

                            # TODO: If this has use_comparison_df = True then we will not have any SOURCE_FILE column
                            # TODO  value we could remove the above lines if use_comparison_df: continue the get and
                            # TODO  print reverse transfer events but later code makes too many assumptions and fails
                            # ...but a result here indicates a series in target, can be extended by source, tag on
                            # extra data
                            # if use_comparison_df:
                            #    new_df['SOURCE_FILE'] = ["TRANSFER_TO_" + series_b_df['SOURCE_FILE'].max()] *
                            # len(series_a_fragids)

                            all_results = all_results.append(new_df, sort=False)
                            # print series_a_df.to_csv()

                    if series_b_df['QUERY_SERIES_SEQID'].tolist()[-1] == 0:

                        result_id += 1
                        series_b_df['RESULT_ID'] = [result_id] * len(series_b_fragids)
                        series_b_df['QUERY_CONTEXT_ID'] = [series_a_df.loc[:, 'CONTEXT_ID'].max()] * \
                            len(series_b_fragids)
                        series_b_df = series_b_df.rename(columns={'SERIES_SEQ_ID': 'RESULT_SERIES_SEQ_ID',
                                                                  'CONTEXT_ID': 'RESULT_CONTEXT_ID',
                                                                  'FRAG_ID': 'RESULT_FRAG_ID',
                                                                  'MOL_ID': 'RESULT_MOL_ID',
                                                                  'ACTIVITY': 'RESULT_MOL_ACTIVITY',
                                                                  'QUERY_SERIES_SEQID': 'QUERY_ORDER'})
                        all_results = all_results.append(series_b_df,  sort=False)
                        # print series_b_df.to_csv()

        # print "Wrapping up with table of shape: ", all_results.shape
        self.logger.info('Auto search FULL done with final table shape %s' % str(all_results.shape))
        return all_results

    def auto_search_write(self, auto_search_result_df, out_csv):
        """
        Writes the result of an auto_search run to chosen filename
        :param auto_search_result_df: pandas df, should come from aut_search method return
        :param out_csv: the file to write results to
        :return: none, writes file to disk
        """
        self.logger.info('Starting auto search and write')
        all_result_ids = auto_search_result_df['RESULT_ID'].unique()

        # validation of df structure
        required_col = ['RESULT_ID', 'SERIES_ID', 'RESULT_SERIES_SEQ_ID', 'QUERY_MOL_ID', 'RESULT_MOL_ID',
                        'RESULT_CONTEXT_ID', 'QUERY_FRAG_ID', 'QUERY_MOL_ID', 'QUERY_CONTEXT_ID', 'RESULT_FRAG_ID',
                        'QUERY_ORDER', 'RESULT_MOL_ACTIVITY']

        for col in required_col:
            if col not in auto_search_result_df.columns:
                raise Exception("Input data table does not have required columns: %s" % col)

        # catch for empty table
        if auto_search_result_df.shape[0] == 0:
            print ("No results found")
            return False

        iteration = 1
        return_df = None

        for result_id in all_result_ids:

            self.logger.info("Result, series ID %s from table size %s: " % (result_id, auto_search_result_df.shape[0]))

            sub_series_df = auto_search_result_df[auto_search_result_df['RESULT_ID'] == result_id]

            # get the original query mol_id_list in it's original query order
            # it can be mis-ordered due to strict_order=False param on the search method
            mol_id_list = list(zip(sub_series_df['QUERY_MOL_ID'].tolist(), sub_series_df['QUERY_ORDER'].tolist()))
            mol_id_list = sorted(mol_id_list, key=lambda xx: xx[1])
            mol_id_list = [x[0] for x in mol_id_list if x[1] > 0]

            self.logger.debug('Merging results to CSV frame for iteration %s and dataframe %s' %
                              (iteration, str(sub_series_df.shape)))

            if iteration == 1:
                return_df = self.return_scored_series_dataframe(mol_id_list, sub_series_df, return_df, append=False)
                self.logger.debug('First iteration, sized at %s' % str(return_df.shape))
                iteration += 1
            else:
                # as above but append=True
                return_df = self.return_scored_series_dataframe(mol_id_list, sub_series_df, return_df, append=True)
                self.logger.debug('Merge operation, sized at %s' % str(return_df.shape))
                iteration += 1

        # return_df = self.enumerate_products(return_df, 'QUERY_MOL_CONTEXT', 'NEW_FRAG_R')

        return_df.to_csv(out_csv, index=False, float_format='%.3f')  # , header=True)
        self.logger.info('Completed write of auto_search results')


#
# unittest everything
#
class _TestMMPSeriesObjectClass(unittest.TestCase):
    """Test class to test the object and methods"""

    @classmethod
    def setUpClass(cls):
        #
        cls.maxDiff = None

        # setup test data
        cls.temp_file_input_csv = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='wt')
        cls.temp_file_input_csv_larger = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='wt')
        cls.temp_file_input_csv_confusion = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='wt')
        cls.temp_file_input_csv_anotherone = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='wt')
        cls.temp_file_input_csv_double = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='wt')
        # output
        cls.temp_file_output_series = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='wt')
        cls.temp_file_output_seriessuggest = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='wt')
        cls.temp_file_output_seriessuggest2 = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='wt')
        cls.temp_file_output_autosearch = tempfile.NamedTemporaryFile(delete=False, encoding='utf-8', mode='wt')
        # yet more for multi series file parse and search
        cls.temp_dir_series_files = tempfile.mkdtemp(suffix='series')

        # setup a logger object
        cls.mmplogger = logging.getLogger('mmpobjectclass_testlogger')
        logging.disable(logging.CRITICAL)

        cls.test_dataset_goldeninput_csv_headers = \
                ['SMILES, ID, PIC50']

        # first 8 smiles are taken from matsy paper JMC 2014
        # https://www.ncbi.nlm.nih.gov/pubmed/24601597
        # https://pubs.acs.org/doi/10.1021/jm500022q
        cls.test_dataset_goldeninput_csv_data = {
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL89089
            'CS(=O)(=O)c1ccc(C2=C(c3ccc(OC)cc3)CC3(C2)CC3)cc1, 001, 7.00': None,
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL88550
            'Fc1c(OC)ccc(c1)C1=C(c2ccc(S(=O)(=O)C)cc2)CC2(C1)CC2, 002, 7.68': None,
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL88899
            'Clc1c(OC)ccc(c1)C1=C(c2ccc(S(=O)(=O)C)cc2)CC2(C1)CC2, 003, 8.51': None,
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL86363
            'Brc1c(OC)ccc(c1)C1=C(c2ccc(S(=O)(=O)C)cc2)CC2(C1)CC2, 004, 8.77': None,
            # non chiral analogues of https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL390742
            'O=C(OC)C1C(c2ccccc2)CC2N(C1CC2)C, 010, 8.30': None,
            # non chiral analogues of https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL222705
            'Fc1ccc(C2CC3N(C(C2C(=O)OC)CC3)C)cc1, 011, 8.00': None,
            # non chiral analogues of https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL375888
            'Clc1ccc(C2CC3N(C(C2C(=O)OC)CC3)C)cc1, 012, 7.77': None,
            # non chiral analogues of https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL97887
            'Brc1ccc(C2CC3N(C(C2C(=O)OC)CC3)C)cc1, 013, 8.89': None
        }

        # a few more I made up to create some identical repeat series from
        # different molecule ID's and therefore trigger frequency data output
        cls.test_dataset_goldeninput_csv_data_extras = {
            # halogen switch analogues of
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL1080278
            'Clc1cc2c([n]c[n]c2Nc2ccccc2)cc1, 081, 6.80': None,
            'Clc1cc2c([n]c[n]c2Nc2ccc(F)cc2)cc1, 082, 7.20': None,
            'Clc1ccc(Nc2[n]c[n]c3c2cc(cc3)Cl)cc1, 083, 7.81': None,
            'Brc1ccc(Nc2[n]c[n]c3c2cc(Cl)cc3)cc1, 084, 8.42': None
        }

        # even more to add confusion:
        # (a) repeat R groups with different activity and different core
        # (b) repeat identical activity mols with different R groups
        # both require more than a simple sort on the series
        cls.test_dataset_goldeninput_csv_data_extra_confusion = {
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL300021/ H analogue, append 00 to id
            'c1ccc(cc1)C2=C(CN(C2=O)c3ccccc3)c4ccc(cc4)S(=O)(=O)C, 30002100, 6.50': None,
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL300021/
            'COc1ccc(cc1)C2=C(CN(C2=O)c3ccccc3)c4ccc(cc4)S(=O)(=O)C, 300021, 6.50': None,
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL298356
            'CS(=O)(=O)c1ccc(cc1)C2=C(C(=O)N(C2)c3ccccc3)c4ccc(F)cc4, 298356, 6.90': None,
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL298356 Cl analogue, append 01 to id
            'CS(=O)(=O)c1ccc(cc1)C2=C(C(=O)N(C2)c3ccccc3)c4ccc(Cl)cc4, 29835601, 7.30': None,
            # https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL298356 Br analogue, append 02 to id
            'CS(=O)(=O)c1ccc(cc1)C2=C(C(=O)N(C2)c3ccccc3)c4ccc(Br)cc4, 29835602, 7.91': None,
            # made up analogue of CHEMBL86363 (004 above) to that should get
            # transferred to any matched series
            'CCCOc1c(OC)ccc(c1)C1=C(c2ccc(S(=O)(=O)C)cc2)CC2(C1)CC2, 005, 8.99': None,
        }

        # this set will only produce a double cut
        # cox inhibitors again but taken from
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3813081/
        # activity pIC50 data is synthetic
        cls.test_dataset_goldeninput_csv_data_extras_doublecut = {
            # compounds from the paper
            'CS(=O)(=O)c1ccc(cc1)c2cc(Br)sc2c3ccc(F)cc3, 697, 8.6': None,
            'CS(=O)(=O)c1ccc(cc1)c2nc([nH]c2c3ccc(F)cc3)C(F)(F)F, 12, 8.7': None,
            'CCSc1nnc(c2ccc(F)cc2)n1c3ccc(cc3)S(=O)(=O)C, 13, 8.8': None,
            'CS(=O)(=O)c1ccc(cc1)C2SCC(=O)N2c3ccc(F)cc3, 14, 8.9': None,
            'Cc1cnc2c(c3ccc(F)cc3)c(nn2c1C)c4ccc(cc4)S(=O)(=O)C, 33, 8.7': None,
            'CS(=O)(=O)c1ccc(cc1)c2nc3CCCCc3n2c4ccc(F)cc4, 36, 8.8': None,
            'CS(=O)(=O)c1ccc(cc1)C2Sc3ccccc3C(=O)N2c4ccc(F)cc4, 35, 8.9': None,
            'CCOc1ccc(cc1)c2c(nn3ncccc23)c4ccc(cc4)S(=O)(=O)C, 32, 9.0': None,
            'COc1ccc(cc1)c2nc3ccccc3nc2c4ccc(cc4)S(=O)(=O)C, 30, 9.1': None,
            # to make a matched series with this we create synthetic data
            # below are copies of above first 4 compounds but id's now prefixed with 90xxx
            # and the activity data incremented by 1.0
            # and the p-Phenyl F is changed to a p-Phenyl Br
            'CS(=O)(=O)c1ccc(cc1)c2cc(Br)sc2c3ccc(Br)cc3, 90697, 9.6': None,
            'CS(=O)(=O)c1ccc(cc1)c2nc([nH]c2c3ccc(Br)cc3)C(F)(F)F, 9012, 9.7': None,
            'CCSc1nnc(c2ccc(Br)cc2)n1c3ccc(cc3)S(=O)(=O)C, 9013, 9.8': None,
            'CS(=O)(=O)c1ccc(cc1)C2SCC(=O)N2c3ccc(Br)cc3, 9014, 9.9': None
        }

        ##################################
        #
        # write test data to temp files
        #
        ##################################

        # csv file - basic test
        cls.temp_file_input_csv.write(', '.join(cls.test_dataset_goldeninput_csv_headers)+"\n")
        for data in list(cls.test_dataset_goldeninput_csv_data.keys()):
            cls.temp_file_input_csv.write(data+"\n")
        cls.temp_file_input_csv.close()

        # extend the basic test data
        cls.temp_file_input_csv_larger.write(', '.join(cls.test_dataset_goldeninput_csv_headers)+"\n")
        for data in list(cls.test_dataset_goldeninput_csv_data.keys()):
            cls.temp_file_input_csv_larger.write(data+"\n")
        for data in list(cls.test_dataset_goldeninput_csv_data_extras.keys()):
            cls.temp_file_input_csv_larger.write(data+"\n")
        cls.temp_file_input_csv_larger.close()

        # full test with corner cases
        cls.temp_file_input_csv_confusion.write(', '.join(cls.test_dataset_goldeninput_csv_headers)+"\n")
        for data in list(cls.test_dataset_goldeninput_csv_data.keys()):
            cls.temp_file_input_csv_confusion.write(data+"\n")
        for data in list(cls.test_dataset_goldeninput_csv_data_extra_confusion.keys()):
            cls.temp_file_input_csv_confusion.write(data+"\n")
        cls.temp_file_input_csv_confusion.close()

        # double
        cls.temp_file_input_csv_double.write(', '.join(cls.test_dataset_goldeninput_csv_headers) + "\n")
        for data in list(cls.test_dataset_goldeninput_csv_data_extras_doublecut.keys()):
            cls.temp_file_input_csv_double.write(data+"\n")
        cls.temp_file_input_csv_double.close()

        # container for results data
        cls.test_dataset_testresults = {}


    @classmethod
    def tearDownClass(cls):
        """Cleanup for end of all tests"""

        os.remove(cls.temp_file_input_csv.name)
        os.remove(cls.temp_file_input_csv_larger.name)
        os.remove(cls.temp_file_input_csv_confusion.name)
        os.remove(cls.temp_file_output_series.name)
        os.remove(cls.temp_file_output_seriessuggest.name)
        os.remove(cls.temp_file_output_seriessuggest2.name)
        os.remove(cls.temp_file_output_autosearch.name)

    def setUp(cls):
        """Setup object for clean reuse in further tests"""

        cls.temp_file_output_series = tempfile.NamedTemporaryFile(delete=True, encoding='utf-8', mode='wt')
        # create empty mmp object each time
        cls.test_mmp_series_object = MMPSeriesObjectClass(cls.mmplogger)

    def tearDown(cls):
        """Tear down object for clean reuse in further tests"""

        # cls.test_mmp_series_object.clean_out_data_seriesobj()
        # reusable data struct
        cls.test_mmp_series_object.clean_out_data_seriesobj()
        cls.test_dataset_testresults.clear()
        # reusable results file
        # os.remove(cls.temp_file_output_series.name)

    def test_iterator_mmp_series_numeric(cls):
        """Test the generation of basic matched series"""

        cls.test_mmp_series_object.setup_mmp_data_for_mms(cls.temp_file_input_csv.name,
                                                          'SMILES', 'ID', 'PIC50',
                                                          3, 0.50001)

        for series_id, series_seq_id, ctx_id, frag_id_L, molid_L, molid_L_act in \
                cls.test_mmp_series_object._iterator_mmp_series_numeric():

            cls.test_dataset_testresults[(series_id, series_seq_id, ctx_id, frag_id_L, molid_L, molid_L_act)] = None

        #print(cls.test_dataset_testresults)
        cls.assertEqual({(1, 1, 8, 7, 1, 7.0): None, (1, 2, 8, 44, 2, 7.68): None, (1, 3, 8, 69, 3, 8.51): None,
                         (1, 4, 8, 91, 4, 8.77): None, (2, 1, 29, 24, 1, 7.0): None, (2, 2, 29, 34, 2, 7.68): None,
                         (2, 3, 29, 60, 3, 8.51): None, (2, 4, 29, 82, 4, 8.77): None, (3, 1, 107, 156, 12, 7.77): None,
                         (3, 2, 107, 134, 11, 8.0): None, (3, 3, 107, 106, 10, 8.3): None,
                         (3, 4, 107, 175, 13, 8.89): None, (4, 1, 122, 60, 12, 7.77): None,
                         (4, 2, 122, 34, 11, 8.0): None, (4, 3, 122, 24, 10, 8.3): None,
                         (4, 4, 122, 82, 13, 8.89): None},
                         cls.test_dataset_testresults)

    def test_iterator_mmp_series_numeric_confusion(cls):
        """Test the generation of matched series with confused input data"""

        cls.test_mmp_series_object.setup_mmp_data_for_mms(cls.temp_file_input_csv_confusion.name,
                                                          'SMILES', 'ID', 'PIC50',
                                                          3, 0.50001)
        for series_id, series_seq_id, ctx_id, frag_id_L, molid_L, molid_L_act in \
                cls.test_mmp_series_object._iterator_mmp_series_numeric():

            cls.test_dataset_testresults[(series_id, series_seq_id, ctx_id, frag_id_L, molid_L, molid_L_act)] = None

        #print(cls.test_dataset_testresults)
        cls.assertEqual({(1, 1, 8, 7, 1, 7.0): None, (1, 2, 8, 44, 2, 7.68): None, (1, 3, 8, 69, 3, 8.51): None,
                         (1, 4, 8, 91, 4, 8.77): None, (1, 5, 8, 285, 5, 8.99): None, (2, 1, 29, 24, 1, 7.0): None,
                         (2, 2, 29, 34, 2, 7.68): None, (2, 3, 29, 60, 3, 8.51): None, (2, 4, 29, 82, 4, 8.77): None,
                         (2, 5, 29, 286, 5, 8.99): None, (3, 1, 107, 156, 12, 7.77): None,
                         (3, 2, 107, 134, 11, 8.0): None, (3, 3, 107, 106, 10, 8.3): None,
                         (3, 4, 107, 175, 13, 8.89): None, (4, 1, 122, 60, 12, 7.77): None,
                         (4, 2, 122, 34, 11, 8.0): None, (4, 3, 122, 24, 10, 8.3): None,
                         (4, 4, 122, 82, 13, 8.89): None, (5, 1, 194, 106, 30002100, 6.5): None,
                         (5, 1, 194, 7, 300021, 6.5): None, (5, 2, 194, 134, 298356, 6.9): None,
                         (5, 3, 194, 156, 29835601, 7.3): None, (5, 4, 194, 175, 29835602, 7.91): None,
                         (6, 1, 208, 24, 30002100, 6.5): None, (6, 1, 208, 9, 300021, 6.5): None,
                         (6, 2, 208, 34, 298356, 6.9): None, (6, 3, 208, 60, 29835601, 7.3): None,
                         (6, 4, 208, 82, 29835602, 7.91): None},
                        cls.test_dataset_testresults)

    def test_iterator_mmp_series_numeric_double(cls):
        """Test the generation of basic matched series"""

        cls.test_mmp_series_object.setup_mmp_data_for_mms(cls.temp_file_input_csv_double.name,
                                                          'SMILES', 'ID', 'PIC50',
                                                          3, 0.45)

        for series_id, series_seq_id, ctx_id, frag_id_L, molid_L, molid_L_act in \
                cls.test_mmp_series_object._iterator_mmp_series_numeric(sgl_or_dbl='double'):

            cls.test_dataset_testresults[(series_id, series_seq_id, ctx_id, frag_id_L, molid_L, molid_L_act)] = None

        #print(cls.test_dataset_testresults)
        cls.assertEqual(
            {(1, 1, 290, 13, 697, 8.6): None, (1, 2, 290, 40, 12, 8.7): None, (1, 2, 290, 101, 33, 8.7): None,
             (1, 3, 290, 61, 13, 8.8): None, (1, 3, 290, 118, 36, 8.8): None, (1, 4, 290, 83, 14, 8.9): None,
             (1, 4, 290, 136, 35, 8.9): None, (2, 1, 24516, 13, 90697, 9.6): None, (2, 2, 24516, 224, 9012, 9.7): None,
             (2, 3, 24516, 61, 9013, 9.8): None, (2, 4, 24516, 83, 9014, 9.9): None},
                        cls.test_dataset_testresults)

    def test_generate_store_mmp_series(cls):
        """Test the generation of matched series with confused input data"""

        cls.test_mmp_series_object.setup_mmp_data_for_mms(cls.temp_file_input_csv_confusion.name,
                                                          'SMILES', 'ID', 'PIC50',
                                                          3, 0.50001)

        cls.test_mmp_series_object.generate_store_mmp_series()
        # print cls.test_mmp_series_object.series_df.to_dict()
        cls.assertEqual(len(cls.test_mmp_series_object.series_df.to_dict()['MOL_ID']), 28)

        cls.test_mmp_series_object.generate_store_mmp_series(sgl_or_dbl_or_both='double')
        # print cls.test_mmp_series_object.series_df.to_dict()
        cls.assertEqual(len(cls.test_mmp_series_object.series_df.to_dict()['MOL_ID']), 17)

        cls.test_mmp_series_object.generate_store_mmp_series(sgl_or_dbl_or_both='both')
        # print cls.test_mmp_series_object.series_df.to_dict()
        cls.assertEqual(len(cls.test_mmp_series_object.series_df.to_dict()['MOL_ID']), 45)

    def test_return_series_matching_idlist_frags(cls):
        """ test query by frag id """

        cls.test_mmp_series_object.setup_mmp_data_for_mms(cls.temp_file_input_csv_confusion.name,
                                                          'SMILES', 'ID', 'PIC50',
                                                          3, 0.50001)

        cls.test_mmp_series_object.generate_store_mmp_series()

        result_ = cls.test_mmp_series_object.return_series_matching_idlist([24, 34, 60],
                                                                           'FRAG_ID')
        #print(result_.to_dict())
        cls.assertDictEqual(
            {'SERIES_ID': {5: 2, 6: 2, 7: 2, 8: 2, 9: 2, 14: 4, 15: 4, 16: 4, 17: 4, 23: 6, 24: 6, 25: 6, 26: 6, 27: 6},
             'SERIES_SEQ_ID': {5: 1, 6: 2, 7: 3, 8: 4, 9: 5, 14: 1, 15: 2, 16: 3, 17: 4, 23: 1, 24: 1, 25: 2, 26: 3,
                               27: 4},
             'CONTEXT_ID': {5: 29, 6: 29, 7: 29, 8: 29, 9: 29, 14: 122, 15: 122, 16: 122, 17: 122, 23: 208, 24: 208,
                            25: 208, 26: 208, 27: 208},
             'FRAG_ID': {5: 24, 6: 34, 7: 60, 8: 82, 9: 286, 14: 60, 15: 34, 16: 24, 17: 82, 23: 24, 24: 9, 25: 34,
                         26: 60, 27: 82},
             'MOL_ID': {5: 1, 6: 2, 7: 3, 8: 4, 9: 5, 14: 12, 15: 11, 16: 10, 17: 13, 23: 30002100, 24: 300021,
                        25: 298356, 26: 29835601, 27: 29835602},
             'ACTIVITY': {5: 7.0, 6: 7.68, 7: 8.51, 8: 8.77, 9: 8.99, 14: 7.77, 15: 8.0, 16: 8.3, 17: 8.89, 23: 6.5,
                          24: 6.5, 25: 6.9, 26: 7.3, 27: 7.91}},
                            result_.to_dict())

    def test_return_series_matching_idlist_frags_strictorder(cls):
        """ test query by frag id - Repeat of above but strict ordering requested """
        cls.test_mmp_series_object.setup_mmp_data_for_mms(cls.temp_file_input_csv_confusion.name,
                                                          'SMILES', 'ID', 'PIC50',
                                                          3, 0.50001)

        cls.test_mmp_series_object.generate_store_mmp_series()

        result_ = cls.test_mmp_series_object.return_series_matching_idlist([24, 34, 60], 'FRAG_ID',
                                                                           strict_ordering=True)
        #print(result_.to_dict())
        cls.assertDictEqual({'SERIES_ID': {5: 2, 6: 2, 7: 2, 8: 2, 9: 2, 23: 6, 24: 6, 25: 6, 26: 6, 27: 6},
                             'SERIES_SEQ_ID': {5: 1, 6: 2, 7: 3, 8: 4, 9: 5, 23: 1, 24: 1, 25: 2, 26: 3, 27: 4},
                             'CONTEXT_ID': {5: 29, 6: 29, 7: 29, 8: 29, 9: 29, 23: 208, 24: 208, 25: 208, 26: 208,
                                            27: 208},
                             'FRAG_ID': {5: 24, 6: 34, 7: 60, 8: 82, 9: 286, 23: 24, 24: 9, 25: 34, 26: 60, 27: 82},
                             'MOL_ID': {5: 1, 6: 2, 7: 3, 8: 4, 9: 5, 23: 30002100, 24: 300021, 25: 298356,
                                        26: 29835601, 27: 29835602},
                             'ACTIVITY': {5: 7.0, 6: 7.68, 7: 8.51, 8: 8.77, 9: 8.99, 23: 6.5, 24: 6.5, 25: 6.9,
                                          26: 7.3, 27: 7.91}},
                            result_.to_dict())

    def test_return_series_matching_idlist_mols(cls):
        """ test search for series by mol id  """

        cls.test_mmp_series_object.setup_mmp_data_for_mms(cls.temp_file_input_csv_confusion.name,
                                                          'SMILES', 'ID', 'PIC50',
                                                          3, 0.50001)

        cls.test_mmp_series_object.generate_store_mmp_series()

        result_ = cls.test_mmp_series_object.return_series_matching_idlist([2, 3, 4], 'MOL_ID')

        #print(result_.to_dict())
        cls.assertDictEqual({'SERIES_ID': {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2},
                             'SERIES_SEQ_ID': {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 1, 6: 2, 7: 3, 8: 4, 9: 5},
                             'CONTEXT_ID': {0: 8, 1: 8, 2: 8, 3: 8, 4: 8, 5: 29, 6: 29, 7: 29, 8: 29, 9: 29},
                             'FRAG_ID': {0: 7, 1: 44, 2: 69, 3: 91, 4: 285, 5: 24, 6: 34, 7: 60, 8: 82, 9: 286},
                             'MOL_ID': {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 1, 6: 2, 7: 3, 8: 4, 9: 5},
                             'ACTIVITY': {0: 7.0, 1: 7.68, 2: 8.51, 3: 8.77, 4: 8.99, 5: 7.0, 6: 7.68, 7: 8.51, 8: 8.77,
                                          9: 8.99}},
                            result_.to_dict())

    def test_search_for_mms_to_extend_molids(cls):
        """"""
        cls.test_mmp_series_object.setup_mmp_data_for_mms(cls.temp_file_input_csv_confusion.name,
                                                          'SMILES', 'ID', 'PIC50',
                                                          3, 0.50001)

        cls.test_mmp_series_object.generate_store_mmp_series()

        # for this dataset, strict_ordering = False will get the same results due to the simplicity of the dataset
        result_ = cls.test_mmp_series_object.search_for_mms_to_extend_molids([30002100, 298356, 29835601, 29835602])
        # NaN becomes nan on return/dict convert so simpler to reduce to zero for comparison
        # result_.fillna(value=0, inplace=True)
        #print(result_.to_dict())
        cls.assertDictEqual(
            {'RESULT_ID': {0: 1, 1: 1, 2: 1, 3: 1, 4: 1}, 'RESULT_SERIES_SEQ_ID': {0: 1, 1: 2, 2: 3, 3: 4, 4: 5},
             'RESULT_MOL_ID': {0: 1, 1: 2, 2: 3, 3: 4, 4: 5}, 'RESULT_CONTEXT_ID': {0: 29, 1: 29, 2: 29, 3: 29, 4: 29},
             'RESULT_FRAG_ID': {0: 24, 1: 34, 2: 60, 3: 82, 4: 286},
             'QUERY_FRAG_ID': {0: 24, 1: 34, 2: 60, 3: 82, 4: 0},
             'QUERY_MOL_ID': {0: 30002100, 1: 298356, 2: 29835601, 3: 29835602, 4: 0},
             'RESULT_MOL_ACTIVITY': {0: 7.0, 1: 7.68, 2: 8.51, 3: 8.77, 4: 8.99},
             'QUERY_CONTEXT_ID': {0: 208, 1: 208, 2: 208, 3: 208, 4: 208}}, result_.to_dict())

    def test_iterator_series_pstats(cls):
        """test it"""
        cls.test_mmp_series_object.setup_mmp_data_for_mms(cls.temp_file_input_csv_confusion.name,
                                                          'SMILES', 'ID', 'PIC50',
                                                          3, 0.50001)

        cls.test_mmp_series_object.generate_store_mmp_series()

        result_df = cls.test_mmp_series_object.return_series_matching_idlist([24, 34, 82], 'FRAG_ID')

        results = {}
        for series, observations, total_series_obs, observed_prob, \
            expected_prob, enrichment, p_value, p_value_corrected, significant \
                in cls.test_mmp_series_object._iterator_series_pstats([24, 34, 82], result_df):

            results[series] = [observations, total_series_obs, observed_prob,
                               expected_prob, enrichment, p_value, p_value_corrected, significant]
            # print to regenerate/print the series data line by line
            # print series, results[series]

        #print(results)
        cls.assertDictEqual({(24, 34, 82): [2, 3, 0.6666666666666666, 0.16666666666666666, 4.0, 0.07407407407407407,
                                            0.37037037037037035, False],
                             (34, 24, 82): [1, 3, 0.3333333333333333, 0.16666666666666666, 2.0, 0.42129629629629617,
                                            2.106481481481481, False]},
                            results)

    def test_auto_search_double(cls):
        """Testing the auto search method"""
        cls.test_mmp_series_object.setup_mmp_data_for_mms(cls.temp_file_input_csv_double.name,
                                                          'SMILES', 'ID', 'PIC50',
                                                          3, 0.50001)

        cls.test_mmp_series_object.generate_store_mmp_series(sgl_or_dbl_or_both='double')

        # print cls.test_mmp_series_object.series_df.to_csv()
        result_df = cls.test_mmp_series_object.auto_search(5, 3, strict_ordering=True)

        #print(result_df.to_dict())
        cls.assertEqual({'RESULT_ID': {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1},
                         'SERIES_ID': {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1},
                         'RESULT_SERIES_SEQ_ID': {0: 1, 1: 2, 2: 2, 3: 3, 4: 3, 5: 4, 6: 4},
                         'RESULT_MOL_ID': {0: 697, 1: 12, 2: 33, 3: 13, 4: 36, 5: 14, 6: 35},
                         'RESULT_CONTEXT_ID': {0: 290, 1: 290, 2: 290, 3: 290, 4: 290, 5: 290, 6: 290},
                         'RESULT_FRAG_ID': {0: 13, 1: 40, 2: 104, 3: 63, 4: 121, 5: 86, 6: 139},
                         'QUERY_CONTEXT_ID': {0: 25185, 1: 25185, 2: 25185, 3: 25185, 4: 25185, 5: 25185, 6: 25185},
                         'QUERY_FRAG_ID': {0: 13, 1: 0, 2: 0, 3: 63, 4: 0, 5: 86, 6: 0},
                         'QUERY_MOL_ID': {0: 90697, 1: 0, 2: 0, 3: 9013, 4: 0, 5: 9014, 6: 0},
                         'QUERY_ORDER': {0: 1, 1: 0, 2: 0, 3: 2, 4: 0, 5: 3, 6: 0},
                         'RESULT_MOL_ACTIVITY': {0: 8.6, 1: 8.7, 2: 8.7, 3: 8.8, 4: 8.8, 5: 8.9, 6: 8.9}},
                         result_df.to_dict())


if __name__ == '__main__':
    unittest.main()
