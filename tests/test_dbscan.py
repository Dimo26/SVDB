import unittest
import numpy

from svdb.export_module import DBSCAN


class TestDBSCAN(unittest.TestCase):

    #test that distant points are not merged
    def test_distant_points(self):
        data = numpy.array([[1,1],[1,101]])
        epsilon=100
        m=2
        result=DBSCAN.cluster(data,epsilon,m)
        assert (result[0] == -1 and result[1] == -1)

    #test that close points are merged
    def test_close_points(self):
        data = numpy.array([[1,1],[1,101]])
        epsilon=200
        m=2
        result=DBSCAN.cluster(data,epsilon,m)
        assert (result[0] == 0 and result[1] == 0)

    #test that small clusters smaller than m are not merged
    def test_small_cluster(self):
        data = numpy.array([[1,1],[1,1],[1,101],[1,101]])
        epsilon=100
        m=3
        result=DBSCAN.cluster(data,epsilon,m)
        assert (result[0] == -1 and result[1] == -1)
