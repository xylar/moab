from pymoab import core
from pymoab import types
from pymoab.scd import ScdInterface
from pymoab.hcoord import HomCoord
from subprocess import call
from driver import test_driver, CHECK_EQ
import numpy as np
import os


def test_scds():
    
    bounds = [0,0,0,10,0,0]
    scd_tst(bounds)
    bounds = [0,0,0,10,10,0]
    scd_tst(bounds)
    bounds = [0,0,0,10,10,10]
    scd_tst(bounds)

def scd_tst(bnds):
    assert len(bnds) == 6

    mb = core.Core()
    scd = ScdInterface(mb)
    boxes = scd.find_boxes()
    assert len(boxes) == 0 

    low = HomCoord(bnds[:3])
    high = HomCoord(bnds[3:])

    scdbox = scd.construct_box(low,high)
    hexes = mb.get_entities_by_type(mb.get_root_set(),types.MBHEX)
    assert 1 == len(scd.find_boxes())
    assert bnds[3]*bnds[4]*bnds[5] == len(hexes)

    check_sequence(scdbox, *bnds)
    evaluate_sequence(scdbox)

    
def check_sequence(box, imin, jmin, kmin, imax, jmax, kmax):

    bmin = box.box_min()
    CHECK_EQ(bmin.i(),imin)
    CHECK_EQ(bmin.j(),jmin)
    CHECK_EQ(bmin.k(),kmin)

    bmax = box.box_max()
    CHECK_EQ(bmax.i(),imax)
    CHECK_EQ(bmax.j(),jmax)
    CHECK_EQ(bmax.k(),kmax)

    bsize = box.box_size()

    CHECK_EQ((bmax-bmin+HomCoord([1,1,1,0])), bsize)
    
def evaluate_sequence(box):

    bmin = box.box_min()
    
    bmax = box.box_max()

    start_vert = box.start_vertex()

    for i in range(bmin[0],bmax[0]):
        for j in range(bmin[1],bmax[1]):
            for k in range(bmin[2],bmax[2]):
                #compute value of start vert
                this_vert = start_vert + (i-bmin[0]) + (j-bmin[1])*(bmax[0]-bmin[0]+1) + (k-bmin[2])*(bmax[1]-bmin[1]+1)*(bmax[0]-bmin[0]+1)                
                temp_vert = box.get_vertex([i,j,k])
                CHECK_EQ(temp_vert,this_vert)

                temp_vert2 = box.get_vertex(HomCoord([i,j,k]))
                CHECK_EQ(temp_vert,this_vert)

                CHECK_EQ(box.get_params(this_vert),[i,j,k])

                CHECK_EQ(box.contains(i,j,k),True)

    start_elem = box.start_element()
    
    for i in range(bmin[0],bmax[0]-1):
        for j in range(bmin[1],bmax[1]-1):
            for k in range(bmin[2],bmax[2]-1):
                #compute value of start elem
                this_elem = start_elem + (i-bmin[0]) + (j-bmin[1])*(bmax[0]-bmin[0]) + (k-bmin[2])*(bmax[1]-bmin[1])*(bmax[0]-bmin[0])                
                temp_elem = box.get_element([i,j,k])

                CHECK_EQ(temp_elem,this_elem)

                temp_elem2 = box.get_element(HomCoord([i,j,k]))
                CHECK_EQ(temp_elem,this_elem)

                CHECK_EQ(box.get_params(this_elem),[i,j,k])

                CHECK_EQ(box.contains(i,j,k),True)


if __name__ == "__main__":
    tests = [test_scds,]
    test_driver(tests)
