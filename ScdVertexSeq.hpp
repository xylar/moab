#ifndef SCDVERTEXSEQ
#define SCDVERTEXSEQ

//
// Class: ScdVertexSeq
//
// Purpose: represent a rectangular vertex block of mesh
//
// A ScdVertex represents a rectangular vertex block of mesh, including both vertices and
// the parametric space used to address those vertices.

#include "EntitySequence.hpp"
#include "HomXform.hpp"

class ScdVertexSeq : public VertexEntitySequence
{

private:

    //! parameter min/max, in homogeneous coords ijkh (extra row for stride eventually)
  HomCoord vertexParams[3];

    //! difference between max and min params plus one (i.e. # VERTICES in
    //! each parametric direction)
  int dIJK[3];
  
    //! difference between max and min params (i.e. # VERTEXS in
    //! each parametric direction)
  int dIJKm1[3];

public:

    //! constructor
  ScdVertexSeq(EntitySequenceManager* seq_mgr,
               const MBEntityHandle start_vertex, 
               const int imin, const int jmin, const int kmin,
               const int imax, const int jmax, const int kmax) ;
  
  virtual ~ScdVertexSeq() {};

    //! get handle of vertex at i, j, k
  MBEntityHandle get_vertex(const int i, const int j, const int k) const;

    //! get handle of vertex at homogeneous coordinates
  MBEntityHandle get_vertex(const HomCoord &coords) const;

    //! get the parameters of a given handle; return MB_FAILURE if vhandle not in this
    //! sequence
  MBErrorCode get_params(const MBEntityHandle vhandle,
                          int &i, int &j, int &k) const;
  
    //! get min params for this vertex
  void min_params(int &i, int &j, int &k) const;
  
    //! get max params for this vertex
  void max_params(int &i, int &j, int &k) const;

    //! get the min params
  const HomCoord &min_params() const;

    //! get the max params
  const HomCoord &max_params() const;
  
    //! get the number of vertices in each direction, inclusive
  void param_extents(int &di, int &dj, int &dk) const;

    //! convenience functions for parameter extents
  int i_min() const {return vertexParams[0].hom_coord()[0];}
  int j_min() const {return vertexParams[0].hom_coord()[1];}
  int k_min() const {return vertexParams[0].hom_coord()[2];}
  int i_max() const {return vertexParams[1].hom_coord()[0];}
  int j_max() const {return vertexParams[1].hom_coord()[1];}
  int k_max() const {return vertexParams[1].hom_coord()[2];}

    //! return whether this vseq's parameter space contains these parameters
  bool contains(const HomCoord &coords) const;
  bool contains(const int i, const int j, const int k) const;
  
  MBEntityHandle get_unused_handle();

};

inline MBEntityHandle ScdVertexSeq::get_vertex(const int i, const int j, 
                                                const int k) const
{
  return mStartEntityHandle + (i-i_min()) + (j-j_min())*dIJK[0] + 
    (k-k_min())*dIJK[0]*dIJK[1];
}

inline MBEntityHandle ScdVertexSeq::get_vertex(const HomCoord &coords) const
{
  return get_vertex(coords.hom_coord()[0], coords.hom_coord()[1], coords.hom_coord()[2]);
}

inline MBErrorCode ScdVertexSeq::get_params(const MBEntityHandle vhandle,
                                             int &i, int &j, int &k) const
{
  if (TYPE_FROM_HANDLE(vhandle) != MBVERTEX) return MB_FAILURE;

  int hdiff = vhandle - mStartEntityHandle;

  k = hdiff / (dIJK[0]*dIJK[1]);
  j = (hdiff - (k*dIJK[0]*dIJK[1])) / dIJK[0];
  i = hdiff % dIJK[0];

  k += vertexParams[0].k();
  j += vertexParams[0].j();
  i += vertexParams[0].i();

  return (vhandle >= mStartEntityHandle &&
          i >= i_min() && i <= i_max() &&
          j >= j_min() && j <= j_max() &&
          k >= k_min() && k <= k_max()) ? MB_SUCCESS : MB_FAILURE;
}
  
  //! get min params for this vertex
inline void ScdVertexSeq::min_params(int &i, int &j, int &k) const
{
  i = i_min();
  j = j_min();
  k = k_min();
}

//! get max params for this vertex
inline void ScdVertexSeq::max_params(int &i, int &j, int &k) const
{
  i = i_max();
  j = j_max();
  k = k_max();
}

inline const HomCoord &ScdVertexSeq::min_params() const 
{
  return vertexParams[0];
}

inline const HomCoord &ScdVertexSeq::max_params() const 
{
  return vertexParams[1];
}

  //! get the number of vertices in each direction, inclusive
inline void ScdVertexSeq::param_extents(int &di, int &dj, int &dk) const
{
  di = dIJK[0];
  dj = dIJK[1];
  dk = dIJK[2];
}

inline bool ScdVertexSeq::contains(const HomCoord &coords) const
{
  return (coords >= vertexParams[0] && coords <= vertexParams[1]) ? true : false;
}

inline bool ScdVertexSeq::contains(const int i, const int j, const int k) const
{
  return contains(HomCoord(i, j, k));
}

#endif
