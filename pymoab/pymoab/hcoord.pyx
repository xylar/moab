"""MOAB Homogenous Coordinate class"""

from pymoab cimport moab
from cython.operator cimport dereference as deref
cimport numpy as np
import numpy as np
from .types import check_error


cdef class HomCoord(object):

    def __cinit__(self, args = None):
        if args == None:
            self.inst = new moab.HomCoord()
        elif 4 == len(args):
            self.inst = new moab.HomCoord(args[0],args[1],args[2],args[3])
        elif 3 == len(args):
            self.inst = new moab.HomCoord(args[0],args[1],args[2])
        else:
            raise Exception

    def __richcmp__(self, y, op):
        if op == 2: # == 
            return self.eq(y)
        else:
            assert False
    
    def __del__(self):
        del self.inst

    def __getitem__(self, int index):
        cdef int val = deref(self.inst)[index]
        return val

    def __add__(HomCoord self, HomCoord a):
        cdef moab.HomCoord sum = deref(self.inst) + deref(a.inst)
        cdef HomCoord h = HomCoord()
        del h.inst
        h.inst = new moab.HomCoord(sum)
        return h

    def __sub__(HomCoord self, HomCoord a):
        cdef moab.HomCoord sum = deref(self.inst) - deref(a.inst)
        cdef HomCoord h = HomCoord()
        del h.inst
        h.inst = new moab.HomCoord(sum)
        return h

    def eq(HomCoord self, HomCoord a):
        return deref(self.inst) == deref(a.inst)
    
    def set(self, i, j, k, h = 1):
        assert type(i) is int
        assert type(j) is int
        assert type(k) is int
        assert type(h) is int
        self.inst.set(<int> i, <int> j, <int> k, <int> h)

    def i(self):
        return self.inst.i()

    def j(self):
        return self.inst.j()

    def k(self):
        return self.inst.k()

    def h(self):
        return self.inst.h()
    
    def length_squared(self):
        return self.inst.length_squared()

    def length(self):
        return self.inst.length()

    def normalize(self):
        self.inst.normalize()
