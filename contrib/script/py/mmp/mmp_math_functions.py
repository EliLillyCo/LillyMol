###################################################################
# Summary: Cantor Pairing function and Reverse Cantor used in MMP code
#
# Notes:
# - Allows us to fold a vector of integer values onto a single integer value reversibly
#   http://stackoverflow.com/questions/919612/mapping-two-integers-to-one-in-a-unique-and-deterministic-way
#   http://mathworld.wolfram.com/PairingFunction.html
# - This allows us to create a dict keyed by a single int value, saving memory over a tuple based key
#   https://guillaume.segu.in/blog/code/487/optimizing-memory-usage-in-python-a-case-study/
#   The alternative is to store a tuple of ints (x, y) which is 8x? more expensive than an int?
# - The function is order specific so (int1, int2) -> int3 and (int2, int1) -> int4 (not int3)
#   and should only be applied to positive integer values
# - 2016 changes give 25% reduction in speed over 1K input SMI:
#   from 0.78s -> 0.71s -> 0.59s removing isinstance then sub-functions
#
#############################################################
import math
import unittest


def cantor(x, y):

    # 2016 uncommented this out to speed things up, int input checked elsewhere
    # catch anything other than positive integer
    # if not isinstance(x, int) or x < 0:
    #    raise Exception("Need positive integers as input")
    # if not isinstance(y, int) or y < 0:
    #    raise Exception("Need positive integers as input")

    return int((x + y) * (x + y + 1) / 2 + y)


# def _w(z):
#    n = math.sqrt(8 * z + 1)
#    return int(math.floor((n - 1) / 2))


# def _t(z):
#    w = _w(z)
#    return (w * w + w) / 2


def inv_cantor(z):

    # 2016 uncommented this out to speed things up, int input checked elsewhere
    # catch anything other than positive integer
    # if not isinstance(z, int) or z < 0:
    #    raise Exception("Need positive integer as input")

    # t = _t(z)
    # w = _w(z)
    # 2016 faster without the use of subfunctions
    n = math.sqrt(8 * z + 1)
    w = int(math.floor((n - 1) / 2))
    t = (w * w + w) / 2
    y = z - t
    x = w - y

    return (int(x), int(y))


class _TestAllFunctions(unittest.TestCase):

    def setUp(self):
        self.test_val_1 = 12345
        self.test_val_2 = 98765
        self.test_val_12 = 6172870370

    def test_cantor(self):
        # input (test_val_1, test_val_2) should become test_val_12
        self.assertEqual(cantor(self.test_val_1, self.test_val_2), self.test_val_12)
#        self.assertRaises(Exception, cantor, (-1,1))
#        self.assertRaises(Exception, cantor, ('string',1))
#        self.assertRaises(Exception, cantor, (10.99,1))

    def test_invcantor(self):
        # input value test_val_12 should return a tuple of test_val_1, test_val_2
        self.assertEqual(inv_cantor(self.test_val_12), (self.test_val_1, self.test_val_2))
#        self.assertRaises(Exception, inv_cantor, -1)
#        self.assertRaises(Exception, inv_cantor, 'string')
#        self.assertRaises(Exception, inv_cantor, 10.99)


if __name__ == '__main__':
    unittest.main()
