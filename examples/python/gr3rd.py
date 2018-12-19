#!/usr/bin/env python
# Create a moab h5m file from a gr3 file
#
# this converter assumes the gr3 file is in latitude/longitude format, 
#  it will be converted to MOAB hdf5 based format (h5m)
#  example file can be downloaded from http://ccrm.vims.edu/yinglong/TMP/
#    hgrid.ll.tri
#  or from http://ftp.mcs.anl.gov/pub/fathom/MeshFiles/hgrid.ll.tri
#
# SCHISM manual downloaded from here:
#  http://ccrm.vims.edu/schismweb/SCHISM_v5.6-Manual.pdf page 52-57

import sys
import os
import math
import numpy as np

from optparse import OptionParser

from pymoab import core
from pymoab import types
from pymoab.rng import Range


print "== Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser()
parser.description = "This script takes a gr3 mesh file and generates a MOAB h5m file."

parser.add_option("-g", "--gr3", dest="gr3File", help="gr3 grid file to input.", default="input.gr3", metavar="FILENAME")
parser.add_option("-m", "--moab", dest="moabFile", help="MOAB grid file name used as output.", default="grid.h5m", metavar="FILENAME")
parser.add_option("-d", "--2d", dest="convert2d", help="convert to 2d, ignore 3rd coordinate", default=False)

for option in parser.option_list:
        if option.default != ("NO", "DEFAULT"):
                option.help += (" " if option.help else "") + "[default: %default]"
                
options, args = parser.parse_args()

if not options.gr3File:
        sys.exit('Error: gr3 input grid file is required.  Specify with -g command line argument.')
if not options.moabFile:
        sys.exit('Error: MOAB output grid file is required.  Specify with -m command line argument.')
        

# open the text file with gr3 format
filein = open(options.gr3File, "r")
mb = core.Core()
# Write to output file
line1=filein.readline()
line1=filein.readline()
radius = 6.3781e6 # radius in meters
data = line1.split()
nume= int(data[0])
numv= int(data[1])
print "number of elements: " , nume, " number of vertices: " , numv
coords = np.zeros(numv*3)
hgts   = np.zeros(numv)
if options.convert2d:
  for i in range(0, numv):
    line=filein.readline()
    data = line.split()
    coords[3*i] = float(data[1])
    coords[3*i+1] = float(data[2])
    coords[3*i+2] = 0 # float(data[3])
    hgts[i] = float(data[3])
else:
  for i in range(0, numv):
    line=filein.readline()
    data = line.split()
    lon = float(data[1]) /180. * math.pi
    lat = float(data[2]) /180. * math.pi
    height =  float(data[3]) # this is in meters? 
    radius_proj = (radius+height)/radius # nrmalized to unit radius
    coords[3*i]   = radius_proj * math.cos(lat) * math.cos(lon)
    coords[3*i+1] = radius_proj * math.cos(lat) * math.sin(lon)
    coords[3*i+2] = radius_proj * math.sin(lat)
    hgts[i] = height 
  
connect=np.zeros((nume, 3), dtype='uint64')
print coords[0], coords[1], coords[2], coords[3]
for i in range(0,nume):
  line=filein.readline()
  data = line.split()
  numverts = int (data[1])
  for j in range(0,3):
    connect[i,j] = int(data[2+j])

mb.create_vertices(coords)
mb.create_elements(types.MBTRI, connect)
  
# set global id for triangles
tris = mb.get_entities_by_type(0,types.MBTRI    , False)
data = range(1,nume+1)
global_id_tag = mb.tag_get_handle("GLOBAL_ID",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
gids = np.array(data)
mb.tag_set_data(global_id_tag, tris, gids)
# set global id for vertices
data = range(1,numv+1)
verts = mb.get_entities_by_type(0,types.MBVERTEX, False)
gids = np.array(data)
mb.tag_set_data(global_id_tag, verts, gids)

# set height tag on vertices, to help a little in depth things
height_tag = mb.tag_get_handle("Height",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
mb.tag_set_data(height_tag, verts, hgts)


try:
        #mb.write_file(options.moabFile)
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
            raise(IOError, "Failed to write MOAB file.")
