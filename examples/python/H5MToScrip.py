#!/usr/bin/env python
# Create a moab h5m file from an SCRIP file
# See for details: http://www.earthsystemmodeling.org/esmf_releases/public/ESMF_5_2_0rp1/ESMF_refdoc/node3.html#SECTION03024000000000000000

import sys, os, math
import netCDF4
import numpy as np

from optparse import OptionParser
from pymoab import core
from pymoab import types
from pymoab.rng import Range

print("== Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser()
parser.description = "This script takes a MOAB h5m file and generates a SCRIP file from it."

parser.add_option("-m", "--moab", dest="moabFile", help="MOAB grid file name used as input.", default="grid.h5m", metavar="FILENAME")
parser.add_option("-s", "--scrip", dest="scripFile", help="SCRIP grid file to output.", default="scrip.nc", metavar="FILENAME")
parser.add_option("-w", "--weights", dest="addWeights", help="Add remapping weights to SCRIP file output.", default=False, action="store_true" )

for option in parser.option_list:
        if option.default != ("NO", "DEFAULT"):
                option.help += (" " if option.help else "") + "[default: %default]"
                
options, args = parser.parse_args()

if not options.moabFile:
        sys.exit('Error: MOAB input grid file is required.  Specify with -m command line argument.')
if not options.scripFile:
        sys.exit('Error: SCRIP output grid file is required.  Specify with -s command line argument.')

try:
    fout = netCDF4.Dataset(options.scripFile, 'w', format='NETCDF4')
    fout.set_auto_mask(False)
    #fout.set_always_mask(False)
except:
    raise(IOError, "Failed to open NetCDF file handle.")

mb = core.Core()
try:
    print("MOAB trying to load filename: %s"%options.moabFile)
    mb.load_file(options.moabFile)
    #mb.load_file("out25.h5m")
except:
    raise(IOError, "Failed to load MOAB file.")

def get_entities_corners(type):
    locelems = mb.get_entities_by_type(0, type, False)
    if type == types.MBTRI:
        corners = 3 if len(locelems) else 0
    elif type == types.MBQUAD:
        corners = 4 if len(locelems) else 0
    else: # MBPOLYGON
        corners=0
        for e in locelems:
            conn = mb.get_connectivity(e)
            ic = len(conn)
            corners = ic if ic > corners else corners

    return corners, len(locelems)

elems = mb.get_entities_by_dimension(0, 2)
verts = mb.get_entities_by_dimension(0, 0)

nvpere = 0
npolys, npe = get_entities_corners(types.MBPOLYGON)
nquads, nqe = get_entities_corners(types.MBQUAD)
ntris, nte = get_entities_corners(types.MBTRI)

nvpere = max(npolys,nquads,ntris)
print("Found %d polygons, %d quadrangles and %d triangles, with a maximum vertex/element number = %d"%(npe,nqe,nte,nvpere))

# Dimensions
grid_size=len(elems)
grid_corners=nvpere
grid_rank=2

grid_size_dim = fout.createDimension("grid_size", grid_size)
grid_corners_dim = fout.createDimension("grid_corners", grid_corners)
grid_rank_dim = fout.createDimension("grid_rank", 2)

grid_center_lon = fout.createVariable("grid_center_lon","f4",("grid_size"),zlib=True)
grid_center_lon.units = "degrees"
grid_center_lat = fout.createVariable("grid_center_lat","f4",("grid_size"),zlib=True)
grid_center_lat.units = "degrees"
grid_corner_lon = fout.createVariable("grid_corner_lon","f4",("grid_size", "grid_corners"),zlib=True)
grid_corner_lon.units = "degrees"
grid_corner_lat = fout.createVariable("grid_corner_lat","f4",("grid_size", "grid_corners"),zlib=True)
grid_corner_lat.units = "degrees"

# Write to output file
global_id_tag = mb.tag_get_handle("GLOBAL_ID",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
#
ReferenceTolerance = 1e-10
for e in range(grid_size):
    conn = mb.get_connectivity(elems[e])
    centroid = mb.get_coords(elems[e])
    grid_center_lon[e] = centroid[0]
    grid_center_lat[e] = centroid[1]
    if (e*25)%grid_size == 0:
        print("Completed %d%% of the conversion"%((e*100)/grid_size))
    #for j in range(grid_corners):
    coords = mb.get_coords(conn)
    for j in range(len(conn)):
        vid=j
        # print("\t conn = %d, coords = [%f, %f]"%(vid,coords[3*vid],coords[3*vid+1]))
        # http://paulbourke.net/geometry/transformationprojection/
        # grid_corner_lon[e, j] = coords[3*vid] 
        # grid_corner_lat[e, j] = coords[3*vid+1]

        if (abs(coords[3*vid+2]) < 1.0 - ReferenceTolerance):
            dLon = math.atan2(coords[3*vid+1],coords[3*vid])
            dLat = math.asin(coords[3*vid+2])

            if (dLon < 0.0):
                dLon += 2.0 * math.pi;
            dLon = dLon / math.pi * 180.0
            dLat = dLat / math.pi * 180.0
        elif coords[3*vid+2] > 0:
            dLon = 0.0
            dLat = 90.0
        else:
            dLon = 0.0
            dLat = -90.0

        grid_corner_lon[e, j] = dLon 
        grid_corner_lat[e, j] = dLat


if options.addWeights:
    # Let us read the weigts from the H5M file as well
    src_grid_rank_dim = fout.createDimension("src_grid_rank", 1)
    dst_grid_rank_dim = fout.createDimension("dst_grid_rank", 1)
    # We actually do not store this info in the h5m file currently
    na_dim = fout.createDimension("n_a", grid_size)
    nb_dim = fout.createDimension("n_b", grid_corners)
    nva_dim = fout.createDimension("nv_a", 2)
    nvb_dim = fout.createDimension("nv_b", 2)

## Let us prepare to write out the file now
fout.meshName = options.scripFile
fout.history = "python H5MToScrip.py -m %s -s %s" % (options.moabFile,options.scripFile)

try:
        fout.close()
        assert os.path.isfile(options.scripFile)
except:
        print("""
        WARNING: .nc file write failed. If netcdf4 and hdf5 support is enabled in this
        build there could be a problem.
        """)
        raise IOError

