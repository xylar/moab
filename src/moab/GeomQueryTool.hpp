#ifndef MOABMC_HPP
#define MOABMC_HPP

#include "MBTagConventions.hpp"
#include "moab/CartVect.hpp"
#include "moab/Range.hpp"
#include "moab/Core.hpp"
#include "moab/GeomUtil.hpp"
#include "moab/FileOptions.hpp"
#include "moab/EntityHandle.hpp"
#include "moab/GeomTopoTool.hpp"

#include <vector>
#include <map>
#include <string>
#include <assert.h>

namespace moab {

class GeomQueryTool
{
public:

  // Constructor
  GeomQueryTool(GeomTopoTool* geomtopotool, bool trace_counting = false,
                double overlap_thickness = 0., double numerical_precision = 0.001);

  // Destructor
  ~GeomQueryTool();

  ErrorCode initialize();

  class RayHistory {

  public:
    /**
     * Clear this entire history-- logically equivalent to creating a new history,
     * but probably more efficient.
     */
    void reset();

    /**
     * Clear the history up to the most recent intersection.  This should be
     * called when a ray changes direction at the site of a surface crossing,
     * a situation that most commonly occurs at a reflecting boundary.
     */
    void reset_to_last_intersection();

    /**
     * Remove the most recent intersection.  This allows a subsequent call
     * along the same ray to return the same intersection.
     */
    void rollback_last_intersection();

    /**
     * @return the number of surface crossings currently represented by this ray history
     */
    int size() const { return prev_facets.size(); }

  private:
    std::vector<EntityHandle> prev_facets;

    friend class GeomQueryTool;

  };

  /**\brief find the next surface crossing from a given point in a given direction
   *
   * This is the primary method of DagMC, enabling ray tracing through a geometry.
   * Given a volume and a ray, it computes the surface ID and distance to the
   * nearest intersection on that volume.  The caller can compute the location of
   * the intersection by adding the distance to the ray.
   *
   * When a series of calls to this function are made along the same ray (e.g. for
   * the purpose of tracking a ray through several volumes), the optional history
   * argument should be given.  The history prevents previously intersected facets
   * from being intersected again.  A single history should be used as long as a
   * ray is proceeding forward without changing direction.  This situation is
   * sometimes referred to as "streaming."
   *
   * If a ray changes direction at an intersection site, the caller should call
   * reset_to_last_intersection() on the history object before the next ray fire.
   *
   * @param volume The volume to fire the ray at.
   * @param ray_start An array of x,y,z coordinates from which to start the ray.
   * @param ray_dir An array of x,y,z coordinates indicating the direction of the ray.
   *                Must be of unit length.
   * @param next_surf Output parameter indicating the next surface intersected by the ray.
   *                If no intersection is found, will be set to 0.
   * @param next_surf_dist Output parameter indicating distance to next_surf.  If next_surf is
   *                0, this value is undefined and should not be used.
   * @param history Optional RayHistory object.  If provided, the facets in the history are
   *                assumed to not intersect with the given ray.  The facet intersected
   *                by this query will also be added to the history.
   * @param dist_limit Optional distance limit.  If provided and > 0, no intersections at a
   *                distance further than this value will be returned.
   * @param ray_orientation Optional ray orientation. If provided determines intersections
   *                along the normal provided, e.g. if -1 allows intersections back along the
   *                the ray direction, Default is 1, i.e. exit intersections
   * @param stats Optional TrvStats object used to measure performance of underlying OBB
   *              ray-firing query.  See OrientedBoxTreeTool.hpp for details.
   *
   */
  ErrorCode ray_fire(const EntityHandle volume,
                     const double ray_start[3], const double ray_dir[3],
                     EntityHandle& next_surf, double& next_surf_dist,
                     RayHistory* history = NULL, double dist_limit = 0,
                     int ray_orientation = 1,
                     OrientedBoxTreeTool::TrvStats* stats = NULL );

  /**\brief Test if a point is inside or outside a volume
   *
   * This method finds the point on the boundary of the volume that is nearest
   * the test point (x,y,z).  If that point is "close" to a surface, a boundary test
   * is performed based on the normal of the surface at that point and the
   * optional ray direction (u,v,w).
   * @param volume The volume to test
   * @param xyz The location to test for volume containment
   * @param result Set to 0 if xyz it outside volume, 1 if inside, and -1 if on boundary.
   * @param Optional direction to use for underlying ray fire query.  Used to ensure
   *        consistent results when a ray direction is known.  If NULL or {0,0,0} is
   *        given, a random direction will be used.
   * @param history Optional RayHistory object to pass to underlying ray fire query.
   *        The history is not modified by this call.
   */
  ErrorCode point_in_volume(const EntityHandle volume,
                            const double xyz[3],
                            int& result,
                            const double* uvw = NULL,
                            const RayHistory* history = NULL );

  /**\brief Robust test if a point is inside or outside a volume using unit sphere area method
   *
   * This test may be more robust that the standard point_in_volume, but is much slower.
   * It does not detect 'on boundary' situations as point_in_volume does.
   * @param volume The volume to test
   * @param xyz The location to test for volume containment
   * @param result Set to 0 if xyz it outside volume, 1 if inside.
   */
  ErrorCode point_in_volume_slow( const EntityHandle volume, const double xyz[3], int& result );


  /** \brief Given a ray starting at a surface of a volume, check whether the ray enters or exits the volume
   *
   * This function is most useful for rays that change directions at a surface crossing.
   * It can be used to check whether a direction change redirects the ray back into the originating
   * volume.
   *
   * @param volume The volume to test
   * @param surface A surface on volume
   * @param xyz A point location on surface
   * @param uvw A (unit) direction vector
   * @param result Set to 1 if ray is entering volume, or 0 if it is leaving
   * @param history Optional ray history object from a previous call to ray_fire.  If present and non-empty,
   *        the history is used to look up the surface facet at which the ray begins.  Absent a
   *        history, the facet nearest to xyz will be looked up.  The history should always be
   *        provided if available, as it avoids the computational expense of a nearest-facet query.
   */
  ErrorCode test_volume_boundary( const EntityHandle volume, const EntityHandle surface,
                                  const double xyz[3], const double uvw[3], int& result,
                                  const RayHistory* history = NULL );

  /**\brief Find the distance to the point on the boundary of the volume closest to the test point
   *
   * @param volume Volume to query
   * @param point Coordinates of test point
   * @param result Set to the minimum distance from point to a surface in volume
   */
  ErrorCode closest_to_location( EntityHandle volume, const double point[3], double& result);

  /** Calculate the volume contained in a 'volume' */
  ErrorCode measure_volume( EntityHandle volume, double& result );

  /** Calculate sum of area of triangles */
  ErrorCode measure_area( EntityHandle surface, double& result );

  /** Get the sense of surfaces wrt a volume.  Sense values are:
   *  {-1 -> reversed, 0 -> both, 1 -> forward}
   */
  ErrorCode surface_sense( EntityHandle volume,
                             int num_surfaces,
                             const EntityHandle* surfaces,
                             int* senses_out );

  /** Get the sense of a single surface wrt a volume.  Sense values are:
   *  {-1 -> reversed, 0 -> both, 1 -> forward}
   */
  ErrorCode surface_sense( EntityHandle volume, EntityHandle surface, int& sense_out );

  /** Get the normal to a given surface at the point on the surface closest to a given point
   *
   * @param surf Surface on which to get normal
   * @param xyz Point on surf
   * @param angle Set to coordinates of surface normal nearest xyz
   * @param history Optional ray history from a previous call to ray_fire().
   *        If present and non-empty, return the normal
   *        of the most recently intersected facet, ignoring xyz.
   */
  ErrorCode get_angle(EntityHandle surf, const double xyz[3], double angle[3],
                      const RayHistory* history = NULL );

  /** Get the volume on the other side of a surface
   *
   * @param A surface to query
   * @param old_volume A volume on one side of surface
   * @param new_volume Output parameter for volume on the other side of surface
   * @return MB_SUCCESS if new_volume was set successfully, error if not.
   */
  ErrorCode next_vol( EntityHandle surface, EntityHandle old_volume,
                      EntityHandle& new_volume );

private:

  /**\brief determine the point membership when the point is effectively on the boundary
   *
   * Called by point_in_volume when the point is with tolerance of the boundary. Compares the
   * ray direction with the surface normal to determine a volume membership.
   */
  ErrorCode boundary_case( EntityHandle volume, int& result,
                             double u, double v, double w,
                             EntityHandle facet,
                             EntityHandle surface);

  /** get the solid angle projected by a facet on a unit sphere around a point
   *  - used by point_in_volume_slow
   */
  ErrorCode poly_solid_angle( EntityHandle face, const CartVect& point, double& area );

  /**\brief State object used in calls to ray_fire()
   *
   * Storage for the "history" of a ray.  This represents the surface facets
   * that the ray is known to have crossed, which cannot be crossed again
   * as long as the ray does not change direction.  It is intended to be used
   * with a series of consecutive calls to ray_fire(), in which a ray passes
   * over potentially many surfaces.
   */

public:

private:

  GeomTopoTool* geomTopoTool;
  Interface* MBI;
  OrientedBoxTreeTool* obbTreeTool;
  bool counting;
  long long int n_pt_in_vol_calls = 0;
  long long int n_ray_fire_calls = 0;
  double overlapThickness, numericalPrecision;
  Tag senseTag;
  EntityHandle impl_compl_handle;

  // temporary storage so functions don't have to reallocate vectors
  // for ray_fire:
  std::vector<double> distList;
  std::vector<EntityHandle> prevFacetList, surfList, facetList;

  std::vector<double>       disList;
  std::vector<EntityHandle> surList, facList;
  std::vector<int>          dirList;

};

}

#endif
