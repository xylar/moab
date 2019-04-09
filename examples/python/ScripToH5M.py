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
parser.description = "This script takes a scrip file and generates a MOAB h5m file."

parser.add_option("-s", "--scrip", dest="scripFile", help="SCRIP grid file to input.", default="scrip.nc", metavar="FILENAME")
parser.add_option("-m", "--moab", dest="moabFile", help="MOAB grid file name used as output.", default="grid.h5m", metavar="FILENAME")
parser.add_option("-r", "--norad", dest="nouseRadians", help="Specify that coordinates for SCRIP are not to be written out in degrees.", default=False, action="store_true")

for option in parser.option_list:
        if option.default != ("NO", "DEFAULT"):
                option.help += (" " if option.help else "") + "[default: %default]"
                
options, args = parser.parse_args()

if not options.scripFile:
        sys.exit('Error: SCRIP input grid file is required.  Specify with -s command line argument.')
if not options.moabFile:
        sys.exit('Error: MOAB output grid file is required.  Specify with -m command line argument.') 

fin = netCDF4.Dataset(options.scripFile, 'r')

mb = core.Core()

# Write to output file
# Dimensions
grid_size=len(fin.dimensions["grid_size"])
grid_corners=len(fin.dimensions["grid_corners"])
grid_rank=len(fin.dimensions["grid_rank"])

grid_corner_lat = fin.variables['grid_corner_lat'][:]
grid_corner_lon = fin.variables['grid_corner_lon'][:]

nverts=grid_size*grid_corners
# create vertices
coords=np.zeros(nverts*3)

connect=np.zeros((grid_size,grid_corners), dtype='uint64')
for e in range(grid_size):
   for j in range(grid_corners):
      dLon = grid_corner_lon[e, j]
      dLat = grid_corner_lat[e, j]
      
      # Check if we need to convert coordinates and transform
      # with an appropriate projection
      if not options.nouseRadians:
        dLon = dLon / 180.0 * math.pi
        dLat = dLat / 180.0 * math.pi

        if dLat > 0.5 * math.pi:
          dLat = 0.5 * math.pi
        if dLat < -0.5 * math.pi:
          dLat = -0.5 * math.pi

        coords[3*(e*grid_corners+j)]   = math.cos(dLon) * math.cos(dLat)
        coords[3*(e*grid_corners+j)+1] = math.sin(dLon) * math.cos(dLat)
        coords[3*(e*grid_corners+j)+2] = math.sin(dLat)
      else:
        # No projection necessary. Copy coordinates
        coords[3*(e*grid_corners+j)]   = dLon
        coords[3*(e*grid_corners+j)+1] = dLat

      connect[e,j] = e*grid_corners+j+1

# Create vertices based on coordinates
mb.create_vertices(coords)
# Set connectivity for elements
mb.create_elements(types.MBPOLYGON, connect)

# set global id for polygons
polys = mb.get_entities_by_type(0,types.MBPOLYGON, False)
data = list(range(1,grid_size+1))
global_id_tag = mb.tag_get_handle("GLOBAL_ID",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
gids = np.array(data)
mb.tag_set_data(global_id_tag,polys,gids)

# set global id for vertices
data = list(range(1,nverts+1))
verts = mb.get_entities_by_type(0,types.MBVERTEX, False)
gids = np.array(data)
mb.tag_set_data(global_id_tag,verts,gids)

try:
        mb.write_file(options.moabFile)
        assert os.path.isfile(options.moabFile)
except:
        try:
            print("""
            WARNING: .h5m file write failed. If hdf5 support is enabled in this
            build there could be a problem.
            """)
            mb.write_file("outfile.vtk")
            assert os.path.isfile("outfile.vtk")
        except:
            raise IOError

