"""MOAB Structured Mesh Interface"""

cimport numpy as np
import numpy as np

from pymoab cimport moab
from cython.operator cimport dereference as deref
from libc.stdlib cimport malloc,free
from .range cimport Range
from .core cimport Core
from .hcoord cimport HomCoord
from .types import check_error
from pymoab import types

cdef class ScdParData(object):
    def __cinit__(self):
        self.inst = new moab.ScdParData()

cdef class ScdInterface(object):
    
    def __cinit__(self, Core c):
        self.interface  = <moab.Interface*> c.inst
        self.inst = new moab.ScdInterface(self.interface,False)

    def __del__(self):
        del self.inst

    def construct_box(self,
                      low,
                      high,
                      coords = None,
                      lperiod = None,
                      par_data = None,
                      bint assn_gids = False,
                      int resolve_shared_ents = -1,
                      exceptions = ()):
        cdef moab.ErrorCode err
        cdef ScdBox scdb = ScdBox() 
        cdef HomCoord hl = low
        cdef HomCoord hh = high
        cdef np.ndarray[np.float64_t, ndim=1] c
        cdef double* coords_ptr = NULL
        cdef int num_c = 0
        if coords is not None:
            c = coords
            num_c = len(coords)
            coords_ptr = <double*> c.data
        err = self.inst.construct_box(deref(hl.inst),
                                      deref(hh.inst),
                                      coords_ptr,
                                      num_c,
                                      scdb.inst,
                                      NULL,
                                      NULL,
                                      assn_gids,
                                      resolve_shared_ents)
        check_error(err,exceptions)
        return scdb
        
    def find_boxes(self, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range rng = Range()
        err = self.inst.find_boxes(deref(rng.inst))
        check_error(err, exceptions)
        return rng
    
cdef class ScdBox(object):

    def __cinit__(self):
        self.inst = <moab.ScdBox*> malloc(sizeof(moab.ScdBox))

    def __del__(self):
        del self.inst

    def box_min(self):
        cdef moab.HomCoord bmin = self.inst.box_min()
        cdef HomCoord h = HomCoord()
        del h.inst
        h.inst = new moab.HomCoord(bmin)
        return h

    def box_max(self):
        cdef moab.HomCoord bmax = self.inst.box_max()
        cdef HomCoord h = HomCoord()
        del h.inst
        h.inst = new moab.HomCoord(bmax)
        return h

    def box_size(self):
        cdef moab.HomCoord bsize = self.inst.box_size()
        cdef HomCoord h = HomCoord()
        del h.inst
        h.inst = new moab.HomCoord(bsize)
        return h
    
    def start_vertex(self):
        cdef moab.EntityHandle startv = self.inst.start_vertex()
        return startv

    def get_vertex(self, args):
        cdef moab.EntityHandle vert
        cdef moab.HomCoord mh
        cdef HomCoord h
        cdef int i = 0, j = 0, k = 0
        if isinstance(args,HomCoord):
            h = args
            mh = deref(h.inst)
            vert = self.inst.get_vertex(mh)
            return vert
        elif 3 == len(args):
            i = args[0]
            j = args[1]
            k = args[2]
            vert = self.inst.get_vertex(i,j,k)
            return vert
        else:
            check_error(types.MB_FAILURE)
        
    def get_params(self, moab.EntityHandle entity, exceptions = ()):
        cdef int i = 0, j = 0, k = 0
        cdef moab.ErrorCode err = self.inst.get_params(entity, i, j, k)
        check_error(err,exceptions)
        return [i,j,k]

    def contains(self, int i, int j, int k):
        return self.inst.contains(i, j, k)
