from pymoab import core
from pymoab import types
from pymoab.scd import ScdInterface
from pymoab.hcoord import HomCoord
from subprocess import call
import numpy as np
import os


def test_scds():
    
    bounds = [0,0,0,10,0,0]
    print bounds    
    scd_tst(bounds)
    bounds = [0,0,0,10,10,0]
    print bounds    
    scd_tst(bounds)
    bounds = [0,0,0,10,10,10]
    print bounds
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
    assert bmin.i() == imin
    assert bmin.j() == jmin
    assert bmin.k() == kmin

    bmax = box.box_max()
    assert bmax.i() == imax
    assert bmax.j() == jmax
    assert bmax.k() == kmax

    bsize = box.box_size()

    assert (bmax-bmin+HomCoord([1,1,1,0])) == (bsize)
    
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
                assert temp_vert == this_vert

                temp_vert2 = box.get_vertex(HomCoord([i,j,k]))
                assert temp_vert == this_vert

                assert box.get_params(this_vert) == [i,j,k]

                assert box.contains(i,j,k)
