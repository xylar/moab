"""Example of pymoab use in translating a structured, uniform yt mesh to MOAB.
yt installation is required for this example to work."""

import os,sys
os.environ['LD_LIBRARY_PATH'] += ":/usr/software/moab/dev/pymoab/lib"
print('Python LD_LIBRARY_PATH : %s' % (os.environ['LD_LIBRARY_PATH']))

from ctypes import *
lib1 = cdll.LoadLibrary('/usr/software/moab/dev/pymoab/lib/libMOAB.so')

import argparse
import numpy as np
from pymoab import core,types,topo_util
from pymoab.tag import Tag

def parse_args(): 

    parser = argparse.ArgumentParser() 

    parser.add_argument('filename', type=str, nargs='?', help="MOAB mesh that needs to be smoothed", default="/home/vijaysm/code/fathom/moab/MeshFiles/unittest/surfrandomtris-4part.h5m")
    #parser.add_argument("filename")
    args = parser.parse_args()

    return args


def laplacian_smooth(mb):
    """Generates a uniform grid in a moab instance (mb) meshset which is representative of the yt grid dataset (ds)"""

    #create a meshset for this grid
    rs = mb.get_root_set()
    ldim = 3
    for idim in [3, 2, 1]:
        if len(mb.get_entities_by_dimension(rs, idim)) > 0:
            ldim = idim
            melems = mb.get_entities_by_dimension(rs, 2)
            num_elems = len(melems)
            break

    mverts = mb.get_entities_by_dimension(rs, 0)
    num_verts = len(mverts)

    print "The " + str(ldim) + "-d mesh contains " + str(num_verts) + " vertices and " + str(num_elems) + " elements!"

    # Laplacian smoother parameters
    maxiter = 10
    iter = 0

    adjs = mb.get_adjacencies(melems[0], 0, True)

    # Now let us set the new "smoothed" coordinate to the vertex
    oldvtxcoords = mb.get_coords(mverts)
    newvtxcoords = np.zeros(3 * num_verts)

    mtu = topo_util.MeshTopoUtil(mb)
    while iter < maxiter:
        ivtx = 0
        print 'Laplacian smoothing iteration: ' + str(iter)
        for vtx in mverts:
            #adjs = mb.get_adjacencies(vtx, 0, False, 1)
            adjs = mtu.get_bridge_adjacencies(vtx, bridge_dim=ldim, to_dim=0)
            # print '\t n(adjacencies): ' + str(len(adjs))
            vtxcoords = mb.get_coords(adjs)
            avgcoords = np.zeros(3)
            i=0
            while i < len(adjs):

                avgcoords[0] += vtxcoords[i*3+0]
                avgcoords[1] += vtxcoords[i*3+1]
                avgcoords[2] += vtxcoords[i*3+2]
                i+=1

            avgcoords /= i
            newvtxcoords[ivtx * 3 + 0] = avgcoords[0]
            newvtxcoords[ivtx * 3 + 1] = avgcoords[1]
            newvtxcoords[ivtx * 3 + 2] = avgcoords[2]
            ivtx += 1

        # Now let us set the new "smoothed" coordinate to the vertex
        mb.set_coords(mverts, newvtxcoords)

        coorderr = np.sum(newvtxcoords - oldvtxcoords)
        print("\t Error = %.6e" % coorderr)
        oldvtxcoords = newvtxcoords
        iter += 1
    
    #return the grid meshset
    return rs

def main():

    args = parse_args()
    filename = args.filename
    print filename

    #establish a moab instance for use
    mb = core.Core()

    print "Loading input dataset " + filename.split("/")[-1] + "..."
    mb.load_file(filename)

    laplacian_smooth(mb)

    #write file
    mb.write_file("test.h5m")
    

if __name__ == "__main__": 
    main()
