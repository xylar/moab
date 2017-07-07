"""MOAB Homogenous Coordinate class"""

from pymoab cimport moab
from cython.operator cimport dereference as deref
cimport numpy as np
import numpy as np
from .types import check_error


cdef class HomCoord(object):

    def __cinit__(self, args = None):
        """
        Constructor
        """
        if args == None:
            self.inst = new moab.HomCoord()
        elif 4 == len(args):
            self.inst = new moab.HomCoord(args[0],args[1],args[2],args[3])
        elif 3 == len(args):
            self.inst = new moab.HomCoord(args[0],args[1],args[2])
        else:
            raise ValueError("""
            Incorrect number of args provided to HomCoord constructor.
            """)

    def __richcmp__(self, y, op):
        """
        Comparator
        
        Only equals operator is currently supported.
        """
        if op == 2: # == 
            return self.eq(y)
        else:
            assert False
    
    def __del__(self):
        """
        Destructor
        """
        del self.inst

    def __getitem__(self, int index):
        """
        Index operator
        """
        cdef int val = deref(self.inst)[index]
        if index < 4:
            return val
        else:
            raise StopIteration

    def __add__(HomCoord self, HomCoord a):
        """
        Addition operator

        Parameters
        ----------
        a : HomCoord to add to this HomCoord instance
        
        Returns
        -------
        A new HomCoord which is the sum of 'self' and 'a'
        """
        cdef moab.HomCoord sum = deref(self.inst) + deref(a.inst)
        cdef HomCoord h = HomCoord()
        del h.inst
        h.inst = new moab.HomCoord(sum)
        return h

    def __sub__(HomCoord self, HomCoord a):
        """
        Subtraction operator

        Parameters
        ----------
        a : HomCoord to subtract to this HomCoord instance
        
        Returns
        -------
        A new HomCoord which is the sum of 'self' and 'a'
        """        
        cdef moab.HomCoord sum = deref(self.inst) - deref(a.inst)
        cdef HomCoord h = HomCoord()
        del h.inst
        h.inst = new moab.HomCoord(sum)
        return h

    def eq(HomCoord self, HomCoord a):
        """
        Equals comparator
        """
        return deref(self.inst) == deref(a.inst)
    
    def set(self, i, j, k, h = 1):
        """
        Sets the coordinates of the HomCoord class
        """
        assert type(i) is int
        assert type(j) is int
        assert type(k) is int
        assert type(h) is int
        self.inst.set(<int> i, <int> j, <int> k, <int> h)

    def i(self):
        """
        Retrieves the first coordinate
        """
        return self.inst.i()

    def j(self):
        """
        Retrieves the second coordinate
        """        
        return self.inst.j()

    def k(self):
        """
        Retrieves the third coordinate
        """                
        return self.inst.k()

    def h(self):
        """
        Retrieves the fourth coordinate
        """                        
        return self.inst.h()
    
    def length_squared(self):
        """
        Returns the length of the HomCoord squared
        """                                
        return self.inst.length_squared()

    def length(self):
        """
        Returns the length of the HomCoord
        """                                        
        return self.inst.length()

    def normalize(self):
        """
        Normalize the coordinates
        """
        self.inst.normalize()

    def __str__(self):
        """
        HomCoord as a string
        """
        prefix = "HomCoord: ["
        suffix = "]"
        outstr = prefix
        for val in self:
            outstr += str(val)
            outstr += ", "
        outstr = outstr[:-2]
        outstr += suffix
        return outstr

    def __repr__(self):
        """
        Representation of class as a string
        """
        return self.__str__()
