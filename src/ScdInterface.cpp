#include "moab/ScdInterface.hpp"
#include "moab/Core.hpp"
#include "SequenceManager.hpp"
#include "EntitySequence.hpp"
#include "StructuredElementSeq.hpp"
#include "VertexSequence.hpp"
#include "ScdVertexData.hpp"
#include "assert.h"
#include <iostream>

#define ERRORR(rval, str) {if (MB_SUCCESS != rval)                    \
      {std::cerr << str; return rval; }}

        
namespace moab 
{
    
ScdInterface::ScdInterface(Core *impl, bool boxes) 
        : mbImpl(impl), 
          searchedBoxes(false),
          boxPeriodicTag(0),
          boxDimsTag(0),
          boxSetTag(0)
{
  if (boxes) find_boxes(scdBoxes);
}

  // Destructor
ScdInterface::~ScdInterface() 
{
  for (Range::iterator rit = scdBoxes.begin(); rit != scdBoxes.end(); rit++)
    delete get_scd_box(*rit);

  if (box_set_tag(false)) 
    mbImpl->tag_delete(box_set_tag());

}

Interface *ScdInterface::impl() const
{
  return mbImpl;
}

ErrorCode ScdInterface::find_boxes(std::vector<ScdBox*> &scd_boxes) 
{
  Range tmp_boxes;
  ErrorCode rval = find_boxes(tmp_boxes);
  if (MB_SUCCESS != rval) return rval;

  for (Range::iterator rit = tmp_boxes.begin(); rit != tmp_boxes.end(); rit++) {
    ScdBox *tmp_box = get_scd_box(*rit);
    if (tmp_box) scd_boxes.push_back(tmp_box);
    else rval = MB_FAILURE;
  }

  return rval;
}

ErrorCode ScdInterface::find_boxes(Range &scd_boxes) 
{
  ErrorCode rval = MB_SUCCESS;
  box_dims_tag();
  if (!searchedBoxes) {
    rval = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &boxDimsTag, NULL, 1, 
                                                scdBoxes, Interface::UNION);
    searchedBoxes = true;
    if (MB_SUCCESS != rval) return rval;
  }

  scd_boxes.merge(scdBoxes);
  return rval;
}

ScdBox *ScdInterface::get_scd_box(EntityHandle eh) 
{
  ScdBox *scd_box = NULL;
  if (!box_set_tag(false)) return scd_box;

  mbImpl->tag_get_data(box_set_tag(), &eh, 1, &scd_box);
  return scd_box;
}

ErrorCode ScdInterface::construct_box(HomCoord low, HomCoord high, double *coords, unsigned int num_coords,
                                      ScdBox *& new_box, bool is_periodic_i, bool is_periodic_j) 
{
    // create a rectangular structured mesh block
  ErrorCode rval;

  HomCoord tmp_size = high - low + HomCoord(1, 1, 1, 0);
  if ((tmp_size[1] && num_coords && (int)num_coords <= tmp_size[0]) ||
      (tmp_size[2] && num_coords && (int)num_coords <= tmp_size[0]*tmp_size[1]))
    return MB_FAILURE;

  rval = create_scd_sequence(low, high, MBVERTEX, 0, new_box);
  ERRORR(rval, "Trouble creating scd vertex sequence.");

  if (num_coords && coords) {
      // set the vertex coordinates
    double *xc, *yc, *zc;
    rval = new_box->get_coordinate_arrays(xc, yc, zc);
    ERRORR(rval, "Couldn't get vertex coordinate arrays.");

    unsigned int i = 0;
    for (int kl = low[2]; kl <= high[2]; kl++) {
      for (int jl = low[1]; jl <= high[1]; jl++) {
        for (int il = low[0]; il <= high[0]; il++) {
          xc[i] = coords[3*i];
          if (new_box->box_size()[1])
            yc[i] = coords[3*i+1];
          if (new_box->box_size()[2])
            zc[i] = coords[3*i+2];
          i++;
        }
      }
    }
  }
  
    // create element sequence
  SequenceManager *seq_mgr = mbImpl->sequence_manager();

  EntitySequence *tmp_seq;
  EntityHandle start_ent;

    // construct the sequence
  EntityType this_tp = MBHEX;
  if (1 >= tmp_size[2]) this_tp = MBQUAD;
  if (1 >= tmp_size[2] && 1 >= tmp_size[1]) this_tp = MBEDGE;
  rval = seq_mgr->create_scd_sequence(low, high, this_tp, 0, start_ent, tmp_seq, is_periodic_i, is_periodic_j);
  ERRORR(rval, "Trouble creating scd element sequence.");

  new_box->elem_seq(tmp_seq);
  new_box->start_element(start_ent);

    // add vertex seq to element seq, forward orientation, unity transform
  rval = new_box->add_vbox(new_box,
                             // p1: imin,jmin
                           low, low, 
                             // p2: imax,jmin
                           low + HomCoord(1, 0, 0),
                           low + HomCoord(1, 0, 0),
                             // p3: imin,jmax
                           low + HomCoord(0, 1, 0),
                           low + HomCoord(0, 1, 0));
  ERRORR(rval, "Error constructing structured element sequence.");

    // add the new hexes to the scd box set; vertices were added in call to create_scd_sequence
  Range tmp_range(new_box->start_element(), new_box->start_element() + new_box->num_elements() - 1);
  rval = mbImpl->add_entities(new_box->box_set(), tmp_range);
  ERRORR(rval, "Couldn't add new hexes to box set.");

  return MB_SUCCESS;
}


ErrorCode ScdInterface::create_scd_sequence(HomCoord low, HomCoord high, EntityType tp,
                                            int starting_id, ScdBox *&new_box,
                                            bool is_periodic_i, bool is_periodic_j)
{
  HomCoord tmp_size = high - low + HomCoord(1, 1, 1, 0);
  if ((tp == MBHEX && 1 >= tmp_size[2]) ||
      (tp == MBQUAD && 1 >= tmp_size[1]))
    return MB_TYPE_OUT_OF_RANGE;

  SequenceManager *seq_mgr = mbImpl->sequence_manager();

  EntitySequence *tmp_seq;
  EntityHandle start_ent, scd_set;

    // construct the sequence
  ErrorCode rval = seq_mgr->create_scd_sequence(low, high, tp, starting_id, start_ent, tmp_seq,
                                                is_periodic_i, is_periodic_j);
  if (MB_SUCCESS != rval) return rval;

    // create the set for this rectangle
  rval = create_box_set(low, high, scd_set);
  if (MB_SUCCESS != rval) return rval;

    // make the ScdBox
  new_box = new ScdBox(this, scd_set, tmp_seq);
  if (!new_box) return MB_FAILURE;

    // set the start vertex/element
  Range new_range;
  if (MBVERTEX == tp) {
    new_range.insert(start_ent, start_ent+new_box->num_vertices()-1);
  }
  else {
    new_range.insert(start_ent, start_ent+new_box->num_elements()-1);
  }
  
    // put the entities in the box set
  rval = mbImpl->add_entities(scd_set, new_range);
  if (MB_SUCCESS != rval) return rval;

    // tag the set with the box
  rval = mbImpl->tag_set_data(box_set_tag(), &scd_set, 1, &new_box);
  if (MB_SUCCESS != rval) return rval;

  return MB_SUCCESS;
}

ErrorCode ScdInterface::create_box_set(const HomCoord low, const HomCoord high,
                                       EntityHandle &scd_set, bool is_periodic_i, bool is_periodic_j) 
{
    // create the set and put the entities in it
  ErrorCode rval = mbImpl->create_meshset(MESHSET_SET, scd_set);
  if (MB_SUCCESS != rval) return rval;
  scdBoxes.insert(scd_set);

    // tag the set with parameter extents
  int boxdims[6];
  for (int i = 0; i < 3; i++) boxdims[i] = low[i];
  for (int i = 0; i < 3; i++) boxdims[3+i] = high[i];
  rval = mbImpl->tag_set_data(box_dims_tag(), &scd_set, 1, boxdims);
  if (MB_SUCCESS != rval) return rval;

  int is_periodic[2] = {is_periodic_i, is_periodic_j};
  rval = mbImpl->tag_set_data(box_periodic_tag(), &scd_set, 1, is_periodic);
  if (MB_SUCCESS != rval) return rval;
  
  return rval;
}

Tag ScdInterface::box_periodic_tag(bool create_if_missing) 
{
  if (boxPeriodicTag || !create_if_missing) return boxPeriodicTag;

  ErrorCode rval = mbImpl->tag_get_handle("BOX_PERIODIC", 2, MB_TYPE_INTEGER, 
                                          boxPeriodicTag, MB_TAG_SPARSE|MB_TAG_CREAT);
  if (MB_SUCCESS != rval) return 0;
  return boxPeriodicTag;
}

Tag ScdInterface::box_dims_tag(bool create_if_missing) 
{
  if (boxDimsTag || !create_if_missing) return boxDimsTag;

  ErrorCode rval = mbImpl->tag_get_handle("BOX_DIMS", 6, MB_TYPE_INTEGER, 
                                          boxDimsTag, MB_TAG_SPARSE|MB_TAG_CREAT);
  if (MB_SUCCESS != rval) return 0;
  return boxDimsTag;
}

Tag ScdInterface::box_set_tag(bool create_if_missing) 
{
  if (boxSetTag || !create_if_missing) return boxSetTag;

  ErrorCode rval = mbImpl->tag_get_handle("__BOX_SET", sizeof(ScdBox*), MB_TYPE_OPAQUE,
                                      boxSetTag, MB_TAG_SPARSE|MB_TAG_CREAT);
  if (MB_SUCCESS != rval) return 0;
  return boxSetTag;
}

ScdBox::ScdBox(ScdInterface *sc_impl, EntityHandle box_set,
               EntitySequence *seq1, EntitySequence *seq2) 
        : scImpl(sc_impl), boxSet(box_set), vertDat(NULL), elemSeq(NULL), startVertex(0), startElem(0)
{
  for (int i = 0; i < 6; i++) boxDims[i] = 0;
  for (int i = 0; i < 2; i++) isPeriodic[i] = false;
  VertexSequence *vseq = dynamic_cast<VertexSequence *>(seq1);
  if (vseq) vertDat = dynamic_cast<ScdVertexData*>(vseq->data());
  if (vertDat) {
      // retrieve the parametric space
    for (int i = 0; i < 3; i++) {
      boxDims[i] = vertDat->min_params()[i];
      boxDims[3+i] = vertDat->max_params()[i];
    }
    startVertex = vertDat->start_handle();
  }

  elemSeq = dynamic_cast<StructuredElementSeq *>(seq2);
  if (!elemSeq)
    elemSeq = dynamic_cast<StructuredElementSeq *>(seq1);

  if (elemSeq) {
    if (vertDat) {
        // check the parametric space to make sure it's consistent
      assert(elemSeq->sdata()->min_params() == HomCoord(boxDims, 3) && 
             (elemSeq->sdata()->max_params() + HomCoord(1, 1, 1)) == HomCoord(boxDims, 3));
  
    } 
    else {
        // get the parametric space from the element sequence
      for (int i = 0; i < 3; i++) {
        boxDims[i] = elemSeq->sdata()->min_params()[i];
        boxDims[3+i] = elemSeq->sdata()->max_params()[i];
      }
    }

    startElem = elemSeq->start_handle();

    elemSeq->is_periodic(isPeriodic);
  }

  assert(vertDat || elemSeq);
  
  boxSize = HomCoord(boxDims+3, 3) - HomCoord(boxDims, 3) + HomCoord(1, 1, 1);
  boxSizeIJ = (boxSize[1] ? boxSize[1] : 1) * boxSize[0];
  boxSizeIM1 = boxSize[0]-(isPeriodic[0] ? 0 : 1);
  boxSizeIJM1 = (boxSize[1] ? (boxSize[1]-(isPeriodic[1] ? 0 : 1)) : 1) * boxSizeIM1;
}

ScdBox::~ScdBox() 
{
    // reset the tag on the set
  ScdBox *tmp_ptr = NULL;
  if (boxSet) scImpl->mbImpl->tag_set_data(scImpl->box_set_tag(), &boxSet, 1, &tmp_ptr);
}

EntityHandle ScdBox::get_vertex_from_seq(int i, int j, int k) const
{
  assert(elemSeq);
  return elemSeq->get_vertex(i, j, k);
}

int ScdBox::box_dimension() const
{
  return (startElem ? scImpl->mbImpl->dimension_from_handle(startElem) : -1);
}

ErrorCode ScdBox::add_vbox(ScdBox *vbox,
                           HomCoord from1, HomCoord to1, 
                           HomCoord from2, HomCoord to2,
                           HomCoord from3, HomCoord to3,
                           bool bb_input,
                           const HomCoord &bb_min,
                           const HomCoord &bb_max)
{
  ScdVertexData *dum_data = dynamic_cast<ScdVertexData*>(vbox->vertDat);
  ErrorCode rval = elemSeq->sdata()->add_vsequence(dum_data, from1, to1, from2, to2, from3, to3,
                                                   bb_input, bb_min, bb_max);
  return rval;
}

bool ScdBox::boundary_complete() const
{
  return elemSeq->boundary_complete();
}

ErrorCode ScdBox::get_coordinate_arrays(double *&xc, double *&yc, double *&zc) 
{
  if (!vertDat) return MB_FAILURE;
  
  xc = reinterpret_cast<double*>(vertDat->get_sequence_data(0));
  yc = reinterpret_cast<double*>(vertDat->get_sequence_data(1));
  zc = reinterpret_cast<double*>(vertDat->get_sequence_data(2));
  return MB_SUCCESS;
}
  
ErrorCode ScdBox::get_coordinate_arrays(const double *&xc, const double *&yc, const double *&zc) const
{
  if (!vertDat) return MB_FAILURE;
  xc = reinterpret_cast<const double*>(vertDat->get_sequence_data(0));
  yc = reinterpret_cast<const double*>(vertDat->get_sequence_data(1));
  zc = reinterpret_cast<const double*>(vertDat->get_sequence_data(2));
  return MB_SUCCESS;
}
  
ErrorCode ScdBox::vert_dat(ScdVertexData *vert_dat)
{
  vertDat = vert_dat;
  return MB_SUCCESS;
}

ErrorCode ScdBox::elem_seq(EntitySequence *elem_seq)
{
  elemSeq = dynamic_cast<StructuredElementSeq*>(elem_seq);
  if (elemSeq) elemSeq->is_periodic(isPeriodic);
  
  if (isPeriodic[0])
    boxSizeIM1 = boxSize[0]-(isPeriodic[0] ? 0 : 1);
  if (isPeriodic[0] || isPeriodic[1])
    boxSizeIJM1 = (boxSize[1] ? (boxSize[1]-(isPeriodic[1] ? 0 : 1)) : 1) * boxSizeIM1;

  return (elemSeq ? MB_SUCCESS : MB_FAILURE);
}  

ErrorCode ScdBox::get_params(EntityHandle ent, HomCoord &ijkd) const 
{
    // check first whether this is an intermediate entity, so we know what to do
  int dimension = box_dimension();
  int this_dim = scImpl->impl()->dimension_from_handle(ent);
  
  if ((0 == this_dim && !vertDat) ||
      (this_dim && this_dim == dimension)) {
    assert(elemSeq);
    return elemSeq->get_params(ent, ijkd[0], ijkd[1], ijkd[2]);
  }

  else if (!this_dim && vertDat)
    return vertDat->get_params(ent, ijkd[0], ijkd[1], ijkd[2]);
    
  else return MB_NOT_IMPLEMENTED;
}
    
    //! Get the entity of specified dimension adjacent to parametric element
    /**
     * \param dim Dimension of adjacent entity being requested
     * \param i Parametric coordinates of cell being evaluated
     * \param j Parametric coordinates of cell being evaluated
     * \param k Parametric coordinates of cell being evaluated
     * \param dir Direction (0, 1, or 2), for getting adjacent edges (2d, 3d) or faces (3d) 
     * \param ent EntityHandle of adjacent entity
     * \param create_if_missing If true, creates the entity if it doesn't already exist
     */
ErrorCode ScdBox::get_adj_edge_or_face(int dim, int i, int j, int k, int dir, EntityHandle &ent,
                                       bool create_if_missing) const 
{
    // describe connectivity of sub-element in static array
    // subconnect[dim-1][dir][numv][ijk] where dimensions are:
    // [dim-1]: dim=1 or 2, so this is 0 or 1
    // [dir]: one of 0..2, for ijk directions in a hex
    // [numv]: number of vertices describing sub entity = 2*dim <= 4
    // [ijk]: 3 values for i, j, k
  int subconnect[2][3][4][3] = {
      {{{0, 0, 0}, {1, 0, 0}, {-1, -1, -1}, {-1, -1, -1}}, // i edge
       {{0, 0, 0}, {0, 1, 0}, {-1, -1, -1}, {-1, -1, -1}}, // j edge
       {{0, 0, 0}, {0, 0, 1}, {-1, -1, -1}, {-1, -1, -1}}}, // k edge

      {{{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}}, // i face
       {{0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {0, 0, 1}}, // j face
       {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}}}; // k face

  if (i < boxDims[0] || i > boxDims[3] || j < boxDims[1] || j > boxDims[4] || k < boxDims[2] || k > boxDims[5])
    return MB_FAILURE;
  
    // get the vertices making up this entity
  EntityHandle verts[4];
  for (int ind = 0; ind < 2*dim; ind++) {
    verts[ind] = get_vertex(i+subconnect[dim-1][dir][ind][0],
                            j+subconnect[dim-1][dir][ind][1],
                            k+subconnect[dim-1][dir][ind][2]);
    if (!verts[ind]) return MB_FAILURE;
  }
  
  Range ents;
  ErrorCode rval = scImpl->impl()->get_adjacencies(verts, 2*dim, dim, false, ents);
  if (MB_SUCCESS != rval) return rval;

  if (ents.size() > 1) return MB_FAILURE;
  
  else if (ents.size() == 1) {
    ent = *ents.begin();
  }
  else if (create_if_missing)
    rval = scImpl->impl()->create_element((1 == dim ? MBEDGE : MBQUAD), verts, 2*dim, ent);
    
  return rval;
}
    
} // namespace moab
