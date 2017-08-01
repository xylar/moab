from pymoab import core
from pymoab import types
from pymoab.hcoord import HomCoord
from math import sqrt
import numpy as np
from driver import test_driver, CHECK_EQ


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
        print("Shouldn't be here. Test fails.")
        raise AssertionError

    h = HomCoord([1,2,3,4])
    CHECK_EQ(h.__str__(), "HomCoord: [1, 2, 3, 4]")    
    CHECK_EQ(h.i(),1)
    CHECK_EQ(h.j(),2)
    CHECK_EQ(h.k(),3)
    CHECK_EQ(h.h(),4)
    CHECK_EQ(h.length_squared(),14)
    CHECK_EQ(h.length(),int(sqrt(14)))
    h.normalize()
    CHECK_EQ(h.length(),1)


    h.set(4,3,2,1)
    CHECK_EQ(h.__str__(), "HomCoord: [4, 3, 2, 1]")        
    CHECK_EQ(h.i(),4)
    CHECK_EQ(h.j(),3)
    CHECK_EQ(h.k(),2)
    CHECK_EQ(h.h(),1)

    # testing for possible bug in iterator
    #these should work
    CHECK_EQ(h[0],4)
    CHECK_EQ(h[1],3)
    CHECK_EQ(h[2],2)
    CHECK_EQ(h[3],1)
    try:
        h[4]
    except:
        pass
    else:
        print("Shouldn't be here. Test fails")
        raise AssertionError
    
if __name__ == "__main__":
    tests = [test_homcoord,]
    test_driver(tests)
