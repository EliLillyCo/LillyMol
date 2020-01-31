from mmp.mmp_stats_functions import diffn  # , _mean_diff_invlog
import random
import cProfile
import pstats

# 100 random seeds generated on command line using random.randint() in loop
seeds = [
    246, 507, 469, 491, 60, 863, 500, 541, 23, 383,
    673, 377, 874, 606, 11, 872, 858, 391, 102, 68,
    82, 880, 891, 588, 293, 201, 343, 67, 40, 142,
    477, 967, 590, 620, 798, 513, 497, 817, 88, 789,
    550, 605, 313, 546, 809, 546, 65, 686, 22, 761,
    481, 211, 907, 37, 87, 140, 259, 59, 7, 68, 395,
    922, 514, 468, 234, 941, 621, 154, 368, 817, 567,
    250, 946, 363, 632, 217, 86, 976, 93, 113, 175,
    737, 428, 75, 479, 40, 139, 848, 881, 859, 969,
    367, 739, 119, 409, 885, 102, 100, 584, 654
]

# now generate 100 lists of random length containing random floats
# this will always give the same results because we fix the seed
random.seed(10)

# generate a list of lists to be used in testing:
mega_list = []


def makelist(random_int):
    """Creates a list of random float values of total list length
    random_int (specified as function param)"""
    num_list = []
    for count in range(random_int):
        num_list.append(random.random())
    return num_list


def get_random_lists(minlen, maxlen, list_of_ints):
    """for every int in the list_of_ints code will return a lists of random floats.
    The method is reproducible, so if the input minlen=1 and maxlen=2 for list_of_ints=[123]
    the result will always be the same.  Used for testing the timing of functions that require
    a lot of lists as input."""
    for seed_ in seeds:
        random.seed(seed_)
        random_int = random.randint(minlen, maxlen)
        # print (random_int)
        elements = makelist(random_int)
        yield elements


def create_mega_dataset():
    """Using a set os seeds (a list of int values) a huge list of random length lists of
     random floats is generated. The lists of random length with have a 50:50 distribution
     of list length <3 versus >=3"""
    # generate some long lists from seeds
    get_random_lists(3, 100, seeds)
    # now generate some short length lists
    get_random_lists(1, 2, seeds)

    list_of_lists = []
    for new_list in get_random_lists(3, 100, seeds):
        list_of_lists.append(new_list)

    for new_list in get_random_lists(1, 2, seeds):
        list_of_lists.append(new_list)

    global mega_list
    for _ in range(100):
        mega_list = mega_list + list_of_lists


def run_this_test():
    """This will iterate across all data in mega_list using each list as input to
    the function of choice"""
    for list_to_test in mega_list:
        # Need to switch this betweem True/False
        # _mean_diff_invlog(list_to_test, inc_low_n_vals=False)
        diffn(list_to_test, 50, inc_low_n_vals=True)

# generate some data
create_mega_dataset()

# profile whatever is specified in run_this_test()
cProfile.run('run_this_test()', '/tmp/restats')

# collect results from file and print the function of interest (hard coded)
p = pstats.Stats('/tmp/restats')
p.print_stats('_mean_diff_invlog')
p.print_stats('diffn')
