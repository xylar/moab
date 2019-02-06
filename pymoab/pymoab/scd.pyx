"""MOAB Structured Mesh Interface"""

cimport numpy as np
import numpy as np

from pymoab cimport moab
from cython.operator cimport dereference as deref
from libc.stdlib cimport malloc,free
from libcpp.vector cimport vector
from .rng cimport Range
from .core cimport Core
from .hcoord cimport HomCoord
from .tag cimport Tag
from .types import check_error, MB_TYPE_DOUBLE, MB_FAILURE
from .types import _DTYPE_CONV, _eh_py_type

cdef void* null = NULL

cdef class ScdParData(object):
    def __cinit__(self):
        self.inst = new moab.ScdParData()

cdef class ScdInterface(object):

    def __cinit__(self, Core c):
        """
        Constructor.

        Requires a moab core, c, to operate on.
        """
        self.interface  = <moab.Interface*> c.inst
        self.inst = new moab.ScdInterface(self.interface,False)

    def __del__(self):
        """
        Destructor.
        """
        del self.inst

    def construct_box(self,
                      low,
                      high,
                      coords = None,
                      bint assn_gids = False,
                      int resolve_shared_ents = -1,
                      lperiodic = None,
                      exceptions = ()):
        """
        Construct a new structured mesh box, including both vertices and
        elements.

        Parameter range for vertex box is [low-high], for elements is
        [low-high). Construct quads by passing in low[2] == high[2], and edges
        by passing in low[1] == high[1] and low[2] == high[2]. The result is
        passed back in a pymoab ScdBox class.

        NOTE:
        The actual mesh is retained in MOAB when the ScdBox is destroyed. To
        actually destroy the mesh, call the destroy_mesh function on the ScdBox
        object first, before destroying it.

        Example
        -------
        mb = pymoab.core.Core()
        scd = pymoab.scd.ScdInterface(mb)

        low = pymoab.HomCoord(0,0,0)
        high = pymoab.HomCoord(10,10,10)

        scdbox = scd.construct_box(low, high)

        Parameters
        ----------
        low : pymoab HomCoord class
            coordinate representing the lower corner of the box in parameter
            space
        high : pymoab HomCoord class
            coordinate representing the upper corner of the box in parameter
            space
        coords : 1-D iterable of floats
            coordinates of the box vertices, interleaved (xyzxyz...); if None,
            the coordinates are set to parametric values
        assign_gids : bool (default False)
            if True, assigns 1-based global ids to vertices using the
            GLOBAL_ID_TAG_NAME tag
        resolve_shared_ents : int (default = -1)
            if != -1, resolves shared entities up to and including dimension
            equal to value
        lperiodic : int[3] (default = None)
            indicates whether or not dimensions of the structured mesh
            are locally periodic
        exceptions : tuple (default is empty tuple)
            contains any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        PyMOAB ScdBox class.

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        ValueError
            if provided coordinates are not of the correct type.
        """
        cdef moab.ErrorCode err
        cdef ScdBox scdb = ScdBox()
        cdef HomCoord hl = low
        cdef HomCoord hh = high
        cdef np.ndarray[np.float64_t, ndim=1] c
        cdef double* coords_ptr = NULL
        cdef int num_c = 0
        cdef int lp[3]
        for i in range(3): lp[i] = 0
        if coords is not None:
            c = np.asarray(coords, dtype = _DTYPE_CONV[MB_TYPE_DOUBLE])
            assert c.ndim == 1
            num_c = len(coords)
            coords_ptr = <double*> c.data
        if lperiodic is not None:
            lp = lperiodic
        err = self.inst.construct_box(deref(hl.inst),
                                      deref(hh.inst),
                                      coords_ptr,
                                      num_c,
                                      scdb.inst,
                                      &(lp[0]),
                                      NULL,
                                      assn_gids,
                                      resolve_shared_ents)
        check_error(err,exceptions)
        return scdb

    def box_set_tag(self, create_if_missing = True):
        """
        Retrieves an internal tag used to link EntitySets and
        ScdBox instances.

        Example
        -------
        structured_box = scd.get_scd_box(structured_box_set)

        Parameters
        ----------
        create_if_missing : bool
            Indicates that this tag should be created if it doesn't exist already.

        Returns
        -------
        A PyMOAB tag.

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        """
        cdef Tag tag = Tag()
        tag.inst = self.inst.box_set_tag(create_if_missing)
        if <void*> tag.inst == null:
            check_error(MB_FAILURE)
        else:
            return tag

    def get_scd_box(self, eh, exceptions = ()):
        """
        Returns all structured mesh blocks in a the PyMOAB core instance
        in a PyMOAB Range.

        Example
        -------
        structured_box = scd.get_scd_box(structured_box_set)

        Parameters
        ----------
        eh : MOAB EntityHAndle
            EntityHandle of the structured box set.

        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        A PyMOAB ScdBox (structured box) object.

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        ValueError
            if an EntityHandle is not of the correct type
        """
        cdef ScdBox struct_box = ScdBox()
        struct_box.inst = self.inst.get_scd_box(<unsigned long> eh)
        if <void*> struct_box.inst == null:
            check_error(MB_FAILURE, exceptions)
        else:
            return struct_box

    def find_boxes(self, exceptions = ()):
        """
        Returns all structured mesh blocks in a the PyMOAB core instance
        in a PyMOAB Range.

        Example
        -------
        scdBoxes = scd.find_boxes()

        Parameters
        ----------
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        A range of EntityHandles representing ScdBox the found instances

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        """
        cdef moab.ErrorCode err
        cdef Range rng = Range()
        err = self.inst.find_boxes(deref(rng.inst))
        check_error(err, exceptions)
        return rng

    def get_boxes(self, exceptions = ()):
        """
        Returns all structured mesh blocks created by the ScdInterface.

        Example
        -------
        scdBoxes = scd.get_boxes()

        Parameters
        ----------
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)

        Returns
        -------
        A list of ScdBox objects.

        Raises
        ------
        MOAB ErrorCode
            if a MOAB error occurs
        """
        cdef moab.ErrorCode err
        cdef vector[moab.ScdBox*] vec_boxes
        err = self.inst.get_boxes(vec_boxes)
        check_error(err, exceptions)
        boxes_out = []
        for i in range(vec_boxes.size()):
            new_box = ScdBox()
            new_box.inst = vec_boxes[i] # replace pointer
            boxes_out.append(new_box)
        return boxes_out

cdef class ScdBox(object):

    def __cinit__(self):
        self.inst = <moab.ScdBox*> malloc(sizeof(moab.ScdBox))

    def __del__(self):
        del self.inst

    def box_min(self):
        """
        Returns the lower (min) corner of the ScdBox as a PyMOAB HomCoord class.
        """
        cdef moab.HomCoord bmin = self.inst.box_min()
        cdef HomCoord h = HomCoord()
        del h.inst
        h.inst = new moab.HomCoord(bmin)
        return h

    def box_max(self):
        """
        Returns the upper (max) corner of the ScdBox as a PyMOAB HomCoord class.
        """
        cdef moab.HomCoord bmax = self.inst.box_max()
        cdef HomCoord h = HomCoord()
        del h.inst
        h.inst = new moab.HomCoord(bmax)
        return h

    def box_size(self):
        """
        Returns the parameter extents for the box as a PyMOAB HomCoord class.
        """
        cdef moab.HomCoord bsize = self.inst.box_size()
        cdef HomCoord h = HomCoord()
        del h.inst
        h.inst = new moab.HomCoord(bsize)
        return h

    def num_vertices(self):
        """
        Returns the number of vertices in the box as an integer.
        """
        return self.inst.num_vertices()

    def num_elements(self):
        """
        Returns the number of elements in the box as an integer.
        """
        return self.inst.num_elements()

    def start_vertex(self):
        """
        Returns the EntityHandle (long int) of the first vertex in the structured mesh box.
        """
        cdef moab.EntityHandle startv = self.inst.start_vertex()
        return _eh_py_type(startv)

    def start_element(self):
        """
        Returns the EntityHandle (long int) of the first element in the structured mesh box.
        """
        cdef moab.EntityHandle startv = self.inst.start_element()
        return _eh_py_type(startv)

    def get_vertex(self, args):
        """
        Returns the vertex handle for parameter values i,j,k. These parameter
        values can be provided either as a list of 3 integers or as a PyMOAB
        HomCoord clas with parameter values for the desired vertex.

        Examples
        --------

        scdbox.get_vertex([5,5,5])

        OR

        vertex_hc = pymoab.HomCoord(5,5,5)
        scdbox.get_vertex(vertex_hc)

        Parameters
        ----------
        args : list of 3 ints or a PyMOAB HomCoord class

        Returns
        -------
        Vertex EntityHandle of the specified vertex.

        Raises
        ------
        MOAB ErrorCode
            if args parameter is not a valid type
        """
        cdef moab.EntityHandle vert
        cdef moab.HomCoord mh
        cdef HomCoord h
        cdef int i = 0, j = 0, k = 0
        if isinstance(args,HomCoord):
            h = args
            mh = deref(h.inst)
            vert = self.inst.get_vertex(mh)
            return _eh_py_type(vert)
        elif 3 == len(args):
            i = args[0]
            j = args[1]
            k = args[2]
            vert = self.inst.get_vertex(i,j,k)
            return _eh_py_type(vert)
        else:
            check_error(MB_FAILURE)

    def get_element(self, args):
        """
        Returns the element handle for parameter values i,j,k. These parameter
        values can be provided either as a list of 3 integers or as a PyMOAB
        HomCoord clas with parameter values for the desired element.

        Examples
        --------

        scdbox.get_element([5,5,5])

        OR

        element_hc = pymoab.HomCoord(5,5,5)
        scdbox.get_element(element_hc)

        Parameters
        ----------
        args : list of 3 ints or a PyMOAB HomCoord class

        Returns
        -------
        Element EntityHandle of the specified element.

        Raises
        ------
        MOAB ErrorCode
            if args parameter is not a valid type
        """
        cdef moab.EntityHandle vert
        cdef moab.HomCoord mh
        cdef HomCoord h
        cdef int i = 0, j = 0, k = 0
        if isinstance(args,HomCoord):
            h = args
            mh = deref(h.inst)
            vert = self.inst.get_element(mh)
            return _eh_py_type(vert)
        elif 3 == len(args):
            i = args[0]
            j = args[1]
            k = args[2]
            vert = self.inst.get_element(i,j,k)
            return _eh_py_type(vert)
        else:
            check_error(MB_FAILURE)

    def get_params(self, moab.EntityHandle entity, exceptions = ()):
        """
        Returns the parametric coordinates of the specified entity.

        Examples
        --------

        ijk = scdbox.get_params(vertex_handle)

        OR

        ijk = scdbox.get_params(hex_handle)

        Parameters
        ----------
        entity : MOAB EntityHandle
           entity handle for which parameter values are returned
        exceptions : tuple (default is empty tuple)
            A tuple containing any error types that should
            be ignored. (see pymoab.types module for more info)
        Returns
        -------
        A list of ints containing the parameter values for the specified entity.

        Raises
        ------
        ValueError
            if the entity parameter is not of the right type
        MOAB ErrorCode
            if a MOAB error occurs
        """
        cdef int i = 0, j = 0, k = 0
        cdef moab.ErrorCode err = self.inst.get_params(entity, i, j, k)
        check_error(err,exceptions)
        return [i,j,k]

    def contains(self, int i, int j, int k):
        """
        Returns True if the specified parameter values are within the ScdBox;
        False if not.
        """
        return self.inst.contains(i, j, k)

    def box_set(self):
        """
        Returns the EntityHandle of the set containing the structured box elements.
        """
        return _eh_py_type(self.inst.box_set())
