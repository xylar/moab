#include "iMesh.h"
#include "MBiMesh.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/IntxMesh/Intx2MeshOnSphere.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"

extern "C" void update_tracer(iMesh_Instance instance,
    iBase_EntitySetHandle imesh_euler_set, int * ierr)
{
  using namespace moab;
  const double radius = 1.;
  const double gtol = 1.e-9;
  const bool debug = false;

  Range ents;
  moab::Interface * mb = MOABI;
  *ierr = 1;

  EntityHandle euler_set = (EntityHandle) imesh_euler_set;

  Intx2MeshOnSphere worker(mb);
  worker.set_radius_source_mesh(radius);
  worker.set_radius_destination_mesh(radius);
  worker.set_error_tolerance(gtol);

  EntityHandle covering_lagr_set;

  ErrorCode rval = mb->create_meshset(MESHSET_SET, covering_lagr_set);
  MB_CHK_SET_ERR_RET(rval, "can't create covering set ");

  // we need to update the correlation tag and remote tuples
  rval = worker.create_departure_mesh_2nd_alg(euler_set, covering_lagr_set);
  MB_CHK_SET_ERR_RET(rval, "can't populate covering set ");

  if (debug) {
    rval = mb->write_file("lagr.h5m", 0, 0, &covering_lagr_set, 1);
    MB_CHK_SET_ERR_RET(rval, "can't write covering set ");
  }

  //
  rval = enforce_convexity(mb, covering_lagr_set);
  MB_CHK_SET_ERR_RET(rval, "can't write covering set ");

  EntityHandle outputSet;
  rval = mb->create_meshset(MESHSET_SET, outputSet);
  MB_CHK_SET_ERR_RET(rval, "can't create output set ");

  rval = worker.intersect_meshes(covering_lagr_set, euler_set, outputSet);
  MB_CHK_SET_ERR_RET(rval, "can't intersect ");

  if (debug) {
    rval = mb->write_file("output.vtk", 0, 0, &outputSet, 1);
    MB_CHK_SET_ERR_RET(rval, "can't write covering set ");
  }

  // tagElem is the average computed at each element, from nodal values
  Tag tagElem = 0;
  std::string tag_name2("TracerAverage");
  rval = mb->tag_get_handle(tag_name2.c_str(), 1, MB_TYPE_DOUBLE, tagElem,
      MB_TAG_DENSE | MB_TAG_CREAT);
  MB_CHK_SET_ERR_RET(rval, "can't get tracer tag ");

  // area of the euler element is fixed, store it; it is used to recompute the averages at each
  // time step
  Tag tagArea = 0;
  std::string tag_name4("Area");
  rval = mb->tag_get_handle(tag_name4.c_str(), 1, MB_TYPE_DOUBLE, tagArea,
      MB_TAG_DENSE | MB_TAG_CREAT);
  MB_CHK_SET_ERR_RET(rval, "can't get area tag");

  rval = worker.update_tracer_data(outputSet, tagElem, tagArea);
  MB_CHK_SET_ERR_RET(rval, "can't update tracer ");

  // everything can be deleted now from intx data; polygons, etc.

  *ierr = 0;
  return;
}


