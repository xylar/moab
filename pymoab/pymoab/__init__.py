from pkg_resources import get_distribution, DistributionNotFound
import os.path

""" 
Name
====

Description
-----------

PyMOAB: A Python interface to Argonne National Lab's Mesh Oriented dAtaBase
(MOAB)

PyMOAB provides a means of interactively interrogating, modifying, and generating MOAB mesh files.

Much of the core functionality of MOAB has been implemented in the core.Core
module. Interaction with this is intended to be largely analagous to interaction
with MOAB via its native C++ API though some modifications have been made to
allow for interaction with the interface using the various native Python
constructs such as lists, tuples, etc.

"""

from pkg_resources import get_distribution, DistributionNotFound
import os.path

try:
    _dist = get_distribution('pymoab')
    # Normalize case for Windows systems
    dist_loc = os.path.normcase(_dist.location)
    here = os.path.normcase(__file__)
    if not here.startswith(os.path.join(dist_loc, 'pymoab')):
        # not installed, but there is another version that *is*
        raise DistributionNotFound
except DistributionNotFound:
    __version__ = 'Please install this project with setup.py'
else:
    __version__ = _dist.version
