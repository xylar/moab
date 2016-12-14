from pymoab import core
from pymoab import types
from pymoab.hcoord import HomCoord
from math import sqrt
import numpy as np

def test_homcoord():


    #try default construction
    h = HomCoord()
    #now with some sample args
    h = HomCoord([1,1,1])
    h = HomCoord([1,1,1,1])

    #now a case that should fail
    try:
        h = HomCoord([1])
    except:
        pass
    else:
        print "Shouldn't be here. Test fails."
        raise AssertionError

    h = HomCoord([1,2,3,4])
    assert(1 == h.i())
    assert(2 == h.j())
    assert(3 == h.k())
    assert(4 == h.h())
    print h.length_squared()
    assert(14 == h.length_squared())
    assert(int(sqrt(14)) == h.length())
    h.normalize()
    assert(1 == h.length())

    h.set(4,3,2,1)
    assert(4 == h.i())
    assert(3 == h.j())
    assert(2 == h.k())
    assert(1 == h.h())
