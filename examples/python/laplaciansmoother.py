"""Example of PyMOAB to apply Laplacian smoothing algorithm to improve mesh quality """

import argparse
import numpy as np
from pymoab import core,types,topo_util,skinner
from pymoab.tag import Tag

def parse_args(): 

    parser = argparse.ArgumentParser() 

    parser.add_argument('filename', type=str, nargs='?', help="MOAB mesh that needs to be smoothed", default="../../MeshFiles/unittest/surfrandomtris-4part.h5m")
    parser.add_argument('maxiter', type=int, nargs='?', help="Maximum number of smoothing iterations to apply", default=10)
    args = parser.parse_args()

    return args


def laplacian_smooth(mb, maxiter):
    """Generates a uniform grid in a moab instance (mb) meshset which is representative of the yt grid dataset (ds)"""

    #create a meshset for this grid
    rs = mb.get_root_set()
    ldim = 3
    for idim in [3, 2, 1]:
        if len(mb.get_entities_by_dimension(rs, idim)) > 0:
            ldim = idim
            melems = mb.get_entities_by_dimension(rs, 2)
            break

    mverts = mb.get_entities_by_dimension(rs, 0)

    print("The " + str(ldim) + "-d mesh contains " + str(mverts.size()) + " vertices and " + str(melems.size()) + " elements!")

    # Compute the skin and fix the boundary vertices - these should not be moved
    sknr = skinner.Skinner(mb)
    skin_verts = sknr.find_skin(rs, melems, True, False)
    print("Found " + str(skin_verts.size()) + " boundary vertices in the mesh")

    # Now let us set the new "smoothed" coordinate to the vertex
    oldvtxcoords = mb.get_coords(mverts)
    # size of coords: (3 * mverts.size())
    newvtxcoords = np.copy(oldvtxcoords)

    mtu = topo_util.MeshTopoUtil(mb)
    iter = 0
    while iter < maxiter:
        ivtx = 0
        print('Laplacian smoothing iteration: ' + str(iter))
        for vtx in mverts:
            if vtx in skin_verts:
                ivtx += 1
                continue

            # Not a fixed node - lets apply our smoothing algorithm
            adjs = mtu.get_bridge_adjacencies(vtx, bridge_dim=ldim, to_dim=0)
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

        coorderr = np.linalg.norm(newvtxcoords - oldvtxcoords)
        print("\t L_2 Error in iterate = %.6e" % coorderr)
        oldvtxcoords = np.copy(newvtxcoords)
        iter += 1
    
    #return the grid meshset
    return rs


def main():

    args = parse_args()
    filename = args.filename
    maxiter = args.maxiter

    #establish a moab instance for use
    mb = core.Core()

    print("Loading input dataset " + filename.split("/")[-1] + "...")
    mb.load_file(filename)

    # Apply Laplacian smoothing while fixing the boundary nodes
    laplacian_smooth(mb, maxiter)

    #write file
    mb.write_file("smoothed_mesh.vtk")
    

if __name__ == "__main__": 
    main()

