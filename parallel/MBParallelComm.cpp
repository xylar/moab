#include "MBInterface.hpp"
#include "MBParallelComm.hpp"
#include "MBWriteUtilIface.hpp"
#include "MBReadUtilIface.hpp"
#include "SequenceManager.hpp"
#include "EntitySequence.hpp"
#include "TagServer.hpp"
#include "MBTagConventions.hpp"
#include "MBSkinner.hpp"
#include "MBParallelConventions.h"
#include "MBCore.hpp"
#include "MBError.hpp"
#include "ElementSequence.hpp"
#include "MBCN.hpp"
#include "RangeMap.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>

#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)
const bool debug = false;
const bool debug_packing = false;

#include <math.h>
#include <assert.h>


extern "C" 
{
#include "minmax.h"
#include "gs.h"
#include "errmem.h"
#include "types.h"
#include "sort.h"
#include "tuple_list.h"
}

#ifdef USE_MPI
#include "mpi.h"
#endif

#define INITIAL_BUFF_SIZE 1024

//#define DEBUG_PACKING 1
#ifdef DEBUG_PACKING
unsigned int __PACK_num = 0, __UNPACK_num = 0, __PACK_count = 0, __UNPACK_count = 0;
std::string __PACK_string, __UNPACK_string;

#define PC(n, m) {\
          if (__PACK_num == (unsigned int)n && __PACK_string == m) __PACK_count++;\
          else {\
            if (__PACK_count > 1) std::cerr << " (" << __PACK_count << "x)";\
            __PACK_count = 1; __PACK_string = m; __PACK_num = n;\
            std::cerr << std::endl << "PACK: " << n << m;\
          }}
#define UPC(n, m) {\
          if (__UNPACK_num == (unsigned int)n && __UNPACK_string == m) __UNPACK_count++;\
          else {\
            if (__UNPACK_count > 1) std::cerr << "(" << __UNPACK_count << "x)";\
            __UNPACK_count = 1; __UNPACK_string = m; __UNPACK_num = n;\
            std::cerr << std::endl << "UNPACK: " << n << m;\
          }}
#else
#define PC(n, m)
#define UPC(n, m)
#endif

#define PACK_INT(buff, int_val) {int tmp_val = int_val; PACK_INTS(buff, &tmp_val, 1);}

#define PACK_INTS(buff, int_val, num) {memcpy(buff, int_val, (num)*sizeof(int)); buff += (num)*sizeof(int); PC(num, " ints");}

#define PACK_DBL(buff, dbl_val, num) {memcpy(buff, dbl_val, (num)*sizeof(double)); buff += (num)*sizeof(double); PC(num, " doubles");}

#define PACK_EH(buff, eh_val, num) {memcpy(buff, eh_val, (num)*sizeof(MBEntityHandle)); buff += (num)*sizeof(MBEntityHandle); PC(num, " handles");}

#define PACK_CHAR_64(buff, char_val) {strcpy((char*)buff, char_val); buff += 64; PC(64, " chars");}

#define PACK_VOID(buff, val, num) {memcpy(buff, val, num); buff += num; PC(num, " void");}

#define PACK_RANGE(buff, rng) {int num_subs = num_subranges(rng); PACK_INTS(buff, &num_subs, 1); PC(num_subs, "-subranged range"); \
          for (MBRange::const_pair_iterator cit = rng.const_pair_begin(); cit != rng.const_pair_end(); cit++) { \
            MBEntityHandle eh = (*cit).first; PACK_EH(buff, &eh, 1); \
            eh = (*cit).second; PACK_EH(buff, &eh, 1);}; }

#define UNPACK_INT(buff, int_val) {UNPACK_INTS(buff, &int_val, 1);}

#define UNPACK_INTS(buff, int_val, num) {memcpy(int_val, buff, (num)*sizeof(int)); buff += (num)*sizeof(int); UPC(num, " ints");}

#define UNPACK_DBL(buff, dbl_val, num) {memcpy(dbl_val, buff, (num)*sizeof(double)); buff += (num)*sizeof(double); UPC(num, " doubles");}

#define UNPACK_EH(buff, eh_val, num) {memcpy(eh_val, buff, (num)*sizeof(MBEntityHandle)); buff += (num)*sizeof(MBEntityHandle); UPC(num, " handles");}

#define UNPACK_CHAR_64(buff, char_val) {strcpy(char_val, (char*)buff); buff += 64; UPC(64, " chars");}

#define UNPACK_VOID(buff, val, num) {memcpy(val, buff, num); buff += num; UPC(num, " void");}

#define UNPACK_RANGE(buff, rng) {int num_subs; UNPACK_INTS(buff, &num_subs, 1); UPC(num_subs, "-subranged range"); MBEntityHandle _eh[2]; \
          for (int i = 0; i < num_subs; i++) { UNPACK_EH(buff, _eh, 2); rng.insert(_eh[0], _eh[1]);}}

#define RR(a) if (MB_SUCCESS != result) {\
          dynamic_cast<MBCore*>(mbImpl)->get_error_handler()->set_last_error(a);\
          return result;}

#define RRA(a) if (MB_SUCCESS != result) {\
      std::string tmp_str; mbImpl->get_last_error(tmp_str);\
      tmp_str.append("\n"); tmp_str.append(a);\
      dynamic_cast<MBCore*>(mbImpl)->get_error_handler()->set_last_error(tmp_str.c_str()); \
      return result;}

enum MBMessageTag {MB_MESG_ANY=MPI_ANY_TAG, 
                   MB_MESG_SIZE,
                   MB_MESG_ENTS,
                   MB_MESG_REMOTE_HANDLES_RANGE,
                   MB_MESG_REMOTE_HANDLES_VECTOR
};

MBParallelComm::MBParallelComm(MBInterface *impl, MPI_Comm comm) 
    : mbImpl(impl), procConfig(comm), sharedpTag(0), sharedpsTag(0),
      sharedhTag(0), sharedhsTag(0), pstatusTag(0), ifaceSetsTag(0),
      partitionTag(0)
{
  myBuffer.resize(INITIAL_BUFF_SIZE);

  tagServer = dynamic_cast<MBCore*>(mbImpl)->tag_server();
  sequenceManager = dynamic_cast<MBCore*>(mbImpl)->sequence_manager();

  int flag = 1;
  int retval = MPI_Initialized(&flag);
  if (MPI_SUCCESS != retval || !flag) {
    int argc = 0;
    char **argv = NULL;
    
      // mpi not initialized yet - initialize here
    retval = MPI_Init(&argc, &argv);
  }
}

MBParallelComm::MBParallelComm(MBInterface *impl,
                               std::vector<unsigned char> &tmp_buff, 
                               MPI_Comm comm) 
    : mbImpl(impl), procConfig(comm), sharedpTag(0), sharedpsTag(0),
      sharedhTag(0), sharedhsTag(0), pstatusTag(0), ifaceSetsTag(0),
      partitionTag(0)
{
  myBuffer.swap(tmp_buff);
  int flag = 1;
  int retval = MPI_Initialized(&flag);
  if (MPI_SUCCESS != retval || !flag) {
    int argc = 0;
    char **argv = NULL;
    
      // mpi not initialized yet - initialize here
    retval = MPI_Init(&argc, &argv);
  }
}

//! assign a global id space, for largest-dimension or all entities (and
//! in either case for vertices too)
MBErrorCode MBParallelComm::assign_global_ids(MBEntityHandle this_set,
                                              const int dimension, 
                                              const int start_id,
                                              const bool largest_dim_only,
                                              const bool parallel) 
{
  MBRange entities[4];
  int local_num_elements[4];
  MBErrorCode result;
  for (int dim = 0; dim <= dimension; dim++) {
    if (dim == 0 || !largest_dim_only || dim == dimension) {
      result = mbImpl->get_entities_by_dimension(this_set, dim, entities[dim]); 
      RRA("Failed to get vertices in assign_global_ids.");
    }

      // need to filter out non-locally-owned entities!!!
    MBRange dum_range;
    for (MBRange::iterator rit = entities[dim].begin(); rit != entities[dim].end(); rit++)
      if (mbImpl->handle_utils().rank_from_handle(*rit) != 
          (unsigned int) mbImpl->proc_rank()) 
        dum_range.insert(*rit);
    entities[dim] = entities[dim].subtract(dum_range);
    
    local_num_elements[dim] = entities[dim].size();
  }
  
    // communicate numbers
  std::vector<int> num_elements(procConfig.proc_size()*4);
#ifdef USE_MPI
  if (procConfig.proc_size() > 1 && parallel) {
    int retval = MPI_Alltoall(local_num_elements, 4, MPI_INT,
                              &num_elements[0], procConfig.proc_size()*4, 
                              MPI_INT, procConfig.proc_comm());
    if (0 != retval) return MB_FAILURE;
  }
  else
#endif
    for (int dim = 0; dim < 4; dim++) num_elements[dim] = local_num_elements[dim];
  
    // my entities start at one greater than total_elems[d]
  int total_elems[4] = {start_id, start_id, start_id, start_id};
  
  for (unsigned int proc = 0; proc < procConfig.proc_rank(); proc++) {
    for (int dim = 0; dim < 4; dim++) total_elems[dim] += num_elements[4*proc + dim];
  }
  
    //assign global ids now
  MBTag gid_tag;
  int zero = 0;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), 
                              MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                              &zero, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) return result;
  
  for (int dim = 0; dim < 4; dim++) {
    if (entities[dim].empty()) continue;
    num_elements.resize(entities[dim].size());
    int i = 0;
    for (MBRange::iterator rit = entities[dim].begin(); rit != entities[dim].end(); rit++)
      num_elements[i++] = total_elems[dim]++;
    
    result = mbImpl->tag_set_data(gid_tag, entities[dim], &num_elements[0]); 
    RRA("Failed to set global id tag in assign_global_ids.");
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::send_entities(const int to_proc,
                                          MBRange &orig_ents,
                                          const bool adjs,
                                          const bool tags,
                                          const bool store_remote_handles,
                                          MBRange &final_ents,
                                          bool wait_all) 
{

    // find a spot in the proc buffers vector
  int index = get_buffers(to_proc);
  
    // pack/send the entities to the destination
  MPI_Request req_vec[2] = {MPI_REQUEST_NULL,  // send req (ownerSBuff)
                             MPI_REQUEST_NULL}; // recv req (remote handles)
  MBErrorCode result = pack_send_entities(to_proc, orig_ents, 
                                          adjs, tags, store_remote_handles, false,
                                          ownerSBuffs[index], ownerRBuffs[index],
                                          req_vec[0], req_vec[1], final_ents);
  RRA("Failed to pack-send entities.");

  MPI_Status statuses[2];

    // if we're storing remote handles, process the return message
  if (store_remote_handles) {
    assert(MPI_REQUEST_NULL != req_vec[1]);
    MPI_Wait(&req_vec[1], &statuses[1]);
    MBRange remote_range;
    unsigned char *buff_ptr = &ownerRBuffs[index][0];
    UNPACK_RANGE(buff_ptr, remote_range);
    result = set_remote_data(final_ents, remote_range, to_proc);
    RRA("Trouble setting remote data range on sent entities.");
  }

  if (wait_all) {
      // wait for all to finish
    int success = MPI_Waitall(2, req_vec, statuses);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed in waitall in send_entities.");
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::recv_entities(const int from_proc,
                                          const bool store_remote_handles,
                                          MBRange &final_ents,
                                          bool wait_all) 
{
  MBErrorCode result;
  
    // find a spot in the proc buffers vector
  int index = get_buffers(from_proc);
  
    // recv/unpack entities from the sender
  MPI_Request req = MPI_REQUEST_NULL;
  MPI_Status status;
  int success = MPI_Recv(&ghostRBuffs[index][0], ghostRBuffs[index].size(), 
                         MPI_UNSIGNED_CHAR, from_proc, 
                         MB_MESG_ANY, procConfig.proc_comm(), 
                         &status);
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("Recv failed in recv_ents.");
  }

    // check type, and if it's just the size, resize buffer and re-call recv
  if (MB_MESG_SIZE == status.MPI_TAG) {
    result = recv_size_buff(from_proc,
                            ghostRBuffs[index],
                            req);
    RRA("Failed to resize recv buffer.");
    MPI_Wait(&req, &status);
  }
  
    // ok, got the actual entities, now unpack
  result = recv_unpack_entities(from_proc, store_remote_handles, false,
                                ghostRBuffs[index], ghostSBuffs[index], 
                                req, final_ents);
  RRA("Failed to recv-unpack message.");

  if (wait_all) {
      // wait for last message to finish
    int success = MPI_Wait(&req, &status);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed in wait in recv_entities.");
    }
  }

  return result;
}

int MBParallelComm::get_buffers(int to_proc) 
{
  int ind = -1;
  std::vector<int>::iterator vit = 
    std::find(buffProcs.begin(), buffProcs.end(), to_proc);
  if (vit == buffProcs.end()) {
    ind = buffProcs.size();
    buffProcs.push_back(to_proc);
    ownerSBuffs[ind].resize(INITIAL_BUFF_SIZE);
    ghostRBuffs[ind].resize(INITIAL_BUFF_SIZE);
      // allocate these other buffs in case we're storing remote handles
    ownerRBuffs[ind].resize(INITIAL_BUFF_SIZE);
    ghostSBuffs[ind].resize(INITIAL_BUFF_SIZE);
      
  }
  else ind = vit - buffProcs.begin();
  assert(ind < MAX_SHARING_PROCS);
  return ind;
}

MBErrorCode MBParallelComm::pack_send_entities(const int to_proc,
                                               MBRange &orig_ents,
                                               const bool adjacencies,
                                               const bool tags,
                                               const bool store_remote_handles,
                                               const bool iface_layer,
                                               std::vector<unsigned char> &send_buff,
                                               std::vector<unsigned char> &recv_buff,
                                               MPI_Request &send_req,
                                               MPI_Request &recv_req,
                                               MBRange &final_ents) 
{
#ifndef USE_MPI
  return MB_FAILURE;
#else

  MBErrorCode result = MB_SUCCESS;
  MBRange whole_range;
  int buff_size;
    
  result = pack_buffer(orig_ents, adjacencies, tags, 
                       store_remote_handles, iface_layer, to_proc,
                       final_ents, send_buff, buff_size); 
  RRA("Failed to pack buffer in pack_send.");

    // now that we know how many entities, post a receive for the 
    // remote range, if requested
  int success;
  if (store_remote_handles && iface_layer) {
    recv_buff.resize(final_ents.size()*sizeof(MBEntityHandle) + sizeof(int));
    
    success = MPI_Irecv(&recv_buff[0], recv_buff.size(), MPI_UNSIGNED_CHAR, to_proc, 
                        MB_MESG_REMOTE_HANDLES_VECTOR, procConfig.proc_comm(), 
                        &recv_req);
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }
  else if (store_remote_handles && !iface_layer) {
    recv_buff.resize(2*num_subranges(final_ents)*sizeof(MBEntityHandle) + sizeof(int));
    
    success = MPI_Irecv(&recv_buff[0], recv_buff.size(), MPI_UNSIGNED_CHAR, to_proc, 
                        MB_MESG_REMOTE_HANDLES_RANGE, procConfig.proc_comm(), 
                        &recv_req);
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }

    // if the message is large, send a first message to tell how large
  if (INITIAL_BUFF_SIZE < buff_size) {
    int tmp_buff_size = -buff_size;
    int success = MPI_Send(&tmp_buff_size, sizeof(int), MPI_UNSIGNED_CHAR, 
                           to_proc, MB_MESG_SIZE, procConfig.proc_comm());
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }
    
    // send the buffer
  success = MPI_Isend(&send_buff[0], buff_size, MPI_UNSIGNED_CHAR, to_proc, 
                      MB_MESG_ENTS, procConfig.proc_comm(), &send_req);
  if (success != MPI_SUCCESS) return MB_FAILURE;

  return result;
#endif
}

MBErrorCode MBParallelComm::recv_size_buff(const int from_proc,
                                           std::vector<unsigned char> &recv_buff,
                                           MPI_Request &recv_req) 
{
    // use the received size to resize buffer, then post another irecv
  recv_buff.resize(-(*((int*)&recv_buff[0])));
  int success = MPI_Irecv(&recv_buff[0], recv_buff.size(), MPI_UNSIGNED_CHAR, from_proc, 
                          MB_MESG_ENTS, procConfig.proc_comm(), &recv_req);
  if (MPI_SUCCESS != success) {
    MBErrorCode result = MB_FAILURE;
    RRA("Failed call to Irecv in recv_size_buff.");
  }
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::recv_unpack_entities(const int from_proc,
                                                 const bool store_remote_handles,
                                                 const bool iface_layer,
                                                 std::vector<unsigned char> &recv_buff,
                                                 std::vector<unsigned char> &send_buff,
                                                 MPI_Request &send_req,
                                                 MBRange &recd_ents)
{
#ifndef USE_MPI
  return MB_FAILURE;
#else

  MBErrorCode result = MB_SUCCESS;
  
    // unpack the buffer
  if (iface_layer) {
    std::vector<MBEntityHandle> recd_ents_tmp;
    unsigned char *buff_ptr = &recv_buff[0];
    result = unpack_iface_entities(buff_ptr, from_proc, recd_ents_tmp);
    RRA("Failed to unpack buffer in recv_unpack.");
    if (store_remote_handles) {
      int recd_size = recd_ents_tmp.size()*sizeof(MBEntityHandle) +
          sizeof(int);
      send_buff.resize(recd_size);
      unsigned char *buff_ptr = &send_buff[0];
      PACK_INT(buff_ptr, recd_ents_tmp.size());
      PACK_EH(buff_ptr, &recd_ents_tmp[0], recd_ents_tmp.size());
      int success = MPI_Isend(&send_buff[0], recd_size, MPI_UNSIGNED_CHAR, from_proc, 
                              MB_MESG_REMOTE_HANDLES_VECTOR, procConfig.proc_comm(), &send_req);
      if (success != MPI_SUCCESS) {
        result = MB_FAILURE;
        RRA("Failed to send handles in recv_unpack.");
      }
    }
    std::copy(recd_ents_tmp.begin(), recd_ents_tmp.end(), mb_range_inserter(recd_ents));
  }
  else {
    result = unpack_buffer(&recv_buff[0], store_remote_handles, from_proc, recd_ents); 
    RRA("Failed to unpack buffer in recv_unpack.");
    if (store_remote_handles) {
      int recd_size = 2*num_subranges(recd_ents)*sizeof(MBEntityHandle) + sizeof(int);
      send_buff.resize(recd_size);
      unsigned char *buff_ptr = &send_buff[0];
      PACK_RANGE(buff_ptr, recd_ents);
      int success = MPI_Isend(&send_buff[0], recd_size, MPI_UNSIGNED_CHAR, from_proc, 
                          MB_MESG_REMOTE_HANDLES_RANGE, procConfig.proc_comm(), &send_req);
      if (success != MPI_SUCCESS) {
        result = MB_FAILURE;
        RRA("Failed to send handles in recv_unpack.");
      }
    }
  }

  return result;
#endif
}

MBErrorCode MBParallelComm::broadcast_entities( const int from_proc,
                                                MBRange &entities,
                                                const bool adjacencies,
                                                const bool tags) 
{
#ifndef USE_MPI
  return MB_FAILURE;
#else
  
  MBErrorCode result = MB_SUCCESS;
  int success;
  MBRange whole_range;
  int buff_size;
  
  std::vector<unsigned char> buff;
  if ((int)procConfig.proc_rank() == from_proc) {
    result = pack_buffer( entities, adjacencies, tags, 
                          false, false, -1,
                          whole_range, buff, buff_size); 
    RRA("Failed to compute buffer size in broadcast_entities.");
  }

  success = MPI_Bcast( &buff_size, 1, MPI_INT, from_proc, procConfig.proc_comm() );
  if (MPI_SUCCESS != success)
    return MB_FAILURE;
  
  if (!buff_size) // no data
    return MB_SUCCESS;

  if ((int)procConfig.proc_rank() != from_proc) 
    buff.resize(buff_size);

  success = MPI_Bcast( &buff[0], buff_size, MPI_UNSIGNED_CHAR, from_proc, procConfig.proc_comm() );
  if (MPI_SUCCESS != success)
    return MB_FAILURE;
  
  if ((int)procConfig.proc_rank() != from_proc) {
    result = unpack_buffer(&buff[0], false, from_proc, entities);
    RRA("Failed to unpack buffer in broadcast_entities.");
  }

  return MB_SUCCESS;
#endif
}

MBErrorCode MBParallelComm::pack_buffer(MBRange &orig_ents, 
                                        const bool adjacencies,
                                        const bool tags,
                                        const bool store_remote_handles,
                                        const bool iface_layer,
                                        const int to_proc,
                                        MBRange &final_ents,
                                        std::vector<unsigned char> &buff,
                                        int &buff_size) 
{
    // pack the buffer with the entity ranges, adjacencies, and tags sections
    // 
    // Note: new entities used in subsequent connectivity lists, sets, or tags, 
    //   are referred to as (MBMAXTYPE + index), where index is into vector 
    //   of new entities, 0-based
    //
    // DATA LAYOUT IN BUFFER:
    // . w/ handles ? (0-no, 1-yes, range, 2-yes, vector)
    // ENTITIES:
    // . for all types:
    //   - ent_type -- end of data indicated by ent_type == MBMAXTYPE
    //   - range(ent type) -- #ents = range.size()
    //   - if (ent_type == MBVERTEX) xxx[#ents], yyy[#ents], zzz[#ents]
    //   - else {
    //     . nodes_per_ent
    //     . connect[#ents * nodes_per_ent]
    //   - }
    //   - if (handles) range/vector of remote handles
    // SETS:
    // . range of set handles -- (#sets = range.size())
    // . options[#sets] (unsigned int)
    // . if (unordered) set range 
    // . else if ordered
    //   - #ents in set
    //   - handles[#ents]
    // . #parents
    // . if (#parents) handles[#parents]
    // . #children
    // . if (#children) handles[#children]
    // ADJACENCIES:
    // (no adjs for now)
    // TAGS:
    // . #tags
    // . for each tag
    //   - tag size
    //   - tag type
    //   - data type
    //   - if (default value)
    //     . 1
    //     . default value
    //   - else
    //     . 0
    //   - name (null-terminated string, in char[64])
    //   - range (size = #ents)
    //   - tag_vals[#ents]

  MBErrorCode result;

  std::vector<MBEntityType> ent_types;
  std::vector<MBRange> all_ranges;
  std::vector<int> verts_per_entity;
  MBRange set_range;
  std::vector<MBRange> set_ranges;
  std::vector<MBTag> all_tags;
  std::vector<MBRange> tag_ranges;
  std::vector<int> set_sizes;
  std::vector<unsigned int> options_vec;

  buff_size = 0;
  MBRange::const_iterator rit;
  
    //=================
    // count first...
    //=================

  unsigned char *buff_ptr = NULL;
  
    // entities
  result = pack_entities(orig_ents, rit, final_ents, buff_ptr,
                         buff_size, true, store_remote_handles, 
                         iface_layer, to_proc, ent_types, all_ranges, 
                         verts_per_entity); 
  RRA("Packing entities (count) failed.");
  
    // sets
  result = pack_sets(orig_ents, rit, final_ents, buff_ptr, buff_size, true,
                     store_remote_handles, to_proc, set_range, set_ranges,
                     set_sizes, options_vec); 
  RRA("Packing sets (count) failed.");
  
    // adjacencies
  if (adjacencies) {
    result = pack_adjacencies(orig_ents, rit, final_ents, buff_ptr, 
                              buff_size, true,
                              store_remote_handles, to_proc);
    RRA("Packing adjs (count) failed.");
  }
    
    // tags
  if (tags) {
    result = pack_tags(orig_ents, rit, final_ents, buff_ptr, 
                       buff_size, true, store_remote_handles, to_proc, 
                       all_tags, tag_ranges);
    RRA("Packing tags (count) failed.");
  }

    //=================
    // now do the real packing
    //=================

  assert(0 <= buff_size);
  buff.resize(buff_size);
  int orig_buff_size = buff_size;
  buff_size = 0;
  buff_ptr = &buff[0];
  
    // entities
  result = pack_entities(orig_ents, rit, final_ents, buff_ptr,
                         buff_size, false, store_remote_handles, 
                         iface_layer, to_proc, ent_types, 
                         all_ranges, verts_per_entity); 
  RRA("Packing entities (real) failed.");
#ifdef DEBUG_PACKING
  std::cerr << "pack_entities buffer space: " << buff_ptr - &buff[0] << " bytes." << std::endl;
  unsigned char *tmp_buff = buff_ptr;
#endif  
  
    // sets
  result = pack_sets(orig_ents, rit, final_ents, buff_ptr, buff_size, false,
                     store_remote_handles, to_proc, set_range, set_ranges,
                     set_sizes, options_vec); 
  RRA("Packing sets (real) failed.");
#ifdef DEBUG_PACKING
  std::cerr << "pack_sets buffer space: " << buff_ptr - tmp_buff << " bytes." << std::endl;
  tmp_buff = buff_ptr;
#endif  
  
    // adjacencies
  if (adjacencies) {
    result = pack_adjacencies(orig_ents, rit, final_ents, buff_ptr, 
                              buff_size, false,
                              store_remote_handles, to_proc);
    RRA("Packing adjs (real) failed.");
  }
    
    // tags
  if (tags) {
    result = pack_tags(orig_ents, rit, final_ents, buff_ptr, 
                       buff_size, false, store_remote_handles, to_proc, all_tags,
                       tag_ranges);
    RRA("Packing tags (real) failed.");
#ifdef DEBUG_PACKING
    std::cerr << "pack_tags buffer space: " << buff_ptr - tmp_buff << " bytes." << std::endl;
    tmp_buff = buff_ptr;
#endif  
  }

    // original buffer size might be larger, because some
    // ranges of local handles pack into more compact ranges of
    // remote handles (because they're expressed using MBMAXTYPE plus
    // index in new entity range, and that list of entities is likely
    // to be compact)
  assert(orig_buff_size >= buff_size);

  return result;
}
 
MBErrorCode MBParallelComm::unpack_buffer(unsigned char *buff_ptr,
                                          const bool store_remote_handles,
                                          int from_proc,
                                          MBRange &final_ents)
{
  if (myBuffer.capacity() == 0) return MB_FAILURE;
  
#ifdef DEBUG_PACKING
    unsigned char *tmp_buff = buff_ptr;
#endif  
  MBErrorCode result = unpack_entities(buff_ptr, store_remote_handles,
                                       from_proc, final_ents);
  RRA("Unpacking entities failed.");
#ifdef DEBUG_PACKING
    std::cerr << "unpack_entities buffer space: " << buff_ptr - tmp_buff << " bytes." << std::endl;
    tmp_buff = buff_ptr;
#endif  
  result = unpack_sets(buff_ptr, final_ents, store_remote_handles, 
                       from_proc);
  RRA("Unpacking sets failed.");
#ifdef DEBUG_PACKING
    std::cerr << "unpack_sets buffer space: " << buff_ptr - tmp_buff << " bytes." << std::endl;
    tmp_buff = buff_ptr;
#endif  
  result = unpack_tags(buff_ptr, final_ents, store_remote_handles,
                       from_proc);
  RRA("Unpacking tags failed.");
#ifdef DEBUG_PACKING
    std::cerr << "unpack_tags buffer space: " << buff_ptr - tmp_buff << " bytes." << std::endl;
    tmp_buff = buff_ptr;
#endif  

#ifdef DEBUG_PACKING
  std::cerr << std::endl;
#endif
  
  return MB_SUCCESS;
}

int MBParallelComm::num_subranges(const MBRange &this_range)
{
    // ok, have all the ranges we'll pack; count the subranges
  int num_sub_ranges = 0;
  for (MBRange::const_pair_iterator pit = this_range.const_pair_begin(); 
       pit != this_range.const_pair_end(); pit++)
    num_sub_ranges++;

  return num_sub_ranges;
}

MBErrorCode MBParallelComm::pack_entities(MBRange &entities,
                                          MBRange::const_iterator &start_rit,
                                          MBRange &whole_range,
                                          unsigned char *&buff_ptr,
                                          int &count,
                                          const bool just_count,
                                          const bool store_remote_handles,
                                          const bool iface_layer,
                                          const int to_proc,
                                          std::vector<MBEntityType> &ent_types, 
                                          std::vector<MBRange> &all_ranges, 
                                          std::vector<int> &verts_per_entity) 
{
    // (. w/ handles ? (0-no, 1-yes, range, 2-yes, vector))
    // ENTITIES:
    // . for all types:
    //   - ent_type -- end of data indicated by ent_type == MBMAXTYPE
    //   - #ents
    //   - if (ent_type == MBVERTEX) xxx[#ents], yyy[#ents], zzz[#ents]
    //   - else {
    //     . nodes_per_ent
    //     . connect[#ents * nodes_per_ent]
    //   - }
    //   - if (handles) range/vector of remote handles

  unsigned char *orig_buff_ptr = buff_ptr;
  MBErrorCode result;
  MBWriteUtilIface *wu = NULL;
  if (!just_count) {
    result = mbImpl->query_interface(std::string("MBWriteUtilIface"), 
                                     reinterpret_cast<void**>(&wu));
    RRA("Couldn't get MBWriteUtilIface.");

  }
  
    // pack vertices
  if (just_count) {
      // don't count size of verts until after all entities are included
    ent_types.push_back(MBVERTEX);
    verts_per_entity.push_back(1);
    all_ranges.push_back(entities.subset_by_type(MBVERTEX));
  }
  else {
    if (!all_ranges[0].empty()) {
      PACK_INT(buff_ptr, ((int) MBVERTEX));
      PACK_INT(buff_ptr, ((int) all_ranges[0].size()));
      int num_verts = all_ranges[0].size();
      std::vector<double*> coords(3);
      for (int i = 0; i < 3; i++)
        coords[i] = reinterpret_cast<double*>(buff_ptr + i * num_verts * sizeof(double));

      assert(NULL != wu);
    
      result = wu->get_node_arrays(3, num_verts, all_ranges[0], 0, 0, coords);
      RRA("Couldn't allocate node space.");
      PC(3*num_verts, " doubles");

      buff_ptr += 3 * num_verts * sizeof(double);

      if (store_remote_handles)
        PACK_RANGE(buff_ptr, all_ranges[0]);

#ifdef DEBUG_PACKING
      std::cerr << "Packed " << all_ranges[0].size() << " vertices." << std::endl;
#endif      
    }
  }

    // place an iterator at the first non-vertex entity
  if (!all_ranges[0].empty()) {
    start_rit = entities.find(*all_ranges[0].rbegin());
    start_rit++;
  }
  else {
    start_rit = entities.begin();
  }
  
  MBRange::const_iterator end_rit = start_rit;
  std::vector<MBRange>::iterator allr_it = all_ranges.begin();
  
    // pack entities
  if (just_count) {    

      // get all ranges of entities that have different #'s of vertices or different types
    while (end_rit != entities.end() && TYPE_FROM_HANDLE(*start_rit) != MBENTITYSET) {

        // get the sequence holding this entity
      EntitySequence *seq;
      ElementSequence *eseq;
      result = sequenceManager->find(*start_rit, seq);
      RRA("Couldn't find entity sequence.");

      if (NULL == seq) return MB_FAILURE;
      eseq = dynamic_cast<ElementSequence*>(seq);

        // if type and nodes per element change, start a new range
      if (eseq->type() != *ent_types.rbegin() || 
          (int) eseq->nodes_per_element() != *verts_per_entity.rbegin()) {
        ent_types.push_back(eseq->type());
        verts_per_entity.push_back(eseq->nodes_per_element());
        all_ranges.push_back(MBRange());
        allr_it++;
      }
    
        // get position in entities list one past end of this sequence
      end_rit = entities.lower_bound(start_rit, entities.end(), eseq->end_handle()+1);

        // put these entities in the last range
      std::copy(start_rit, end_rit, mb_range_inserter(*all_ranges.rbegin()));

        // remove any which are already shared
      if (store_remote_handles && !iface_layer) {
        result = remove_nonowned_shared(*all_ranges.rbegin(), to_proc, false);
        RRA("Failed nonowned_shared test.");
      }

      whole_range.merge(*all_ranges.rbegin());
      
        // now start where we last left off
      start_rit = end_rit;
    }

      // update vertex range to cover vertices included in other entities
      // and remove entities already there
    result = mbImpl->get_adjacencies(whole_range, 0, false, all_ranges[0], 
                                     MBInterface::UNION);
    RRA("Failed get_adjacencies.");

    if (store_remote_handles) {
      whole_range = whole_range.subtract(all_ranges[0]);
      result = remove_nonowned_shared(all_ranges[0], to_proc, false);
      RRA("Failed nonowned_shared test.");
    }
    
      // count those data, now that we know which 
      // entities get communicated
    whole_range.merge(all_ranges[0]);
    count += 3 * sizeof(double) * all_ranges[0].size();
    
      // space for the ranges
    std::vector<MBRange>::iterator vit = all_ranges.begin();
    std::vector<int>::iterator iit = verts_per_entity.begin();
    std::vector<MBEntityType>::iterator eit = ent_types.begin();
    for (; vit != all_ranges.end(); vit++, iit++, eit++) {
        // entity type of this range, but passed as int
      count += sizeof(int);
      
        // number of entities
      count += sizeof(int);

        // nodes per entity
      if (iit != verts_per_entity.begin()) count += sizeof(int);

        // connectivity of subrange
      if (iit != verts_per_entity.begin()) {
        count += *iit * sizeof(MBEntityHandle)*(*vit).size();
      }
      
        // remote handles, if desired
      if (store_remote_handles)
        count += sizeof(int) + 2*sizeof(MBEntityHandle)*num_subranges(*vit);
    }

      // extra entity type at end, passed as int
    count += sizeof(int);
  }
  else {
      // for each range beyond the first
    allr_it++;
    std::vector<int>::iterator nv_it = verts_per_entity.begin();
    std::vector<MBEntityType>::iterator et_it = ent_types.begin();
    nv_it++; et_it++;
    
    for (; allr_it != all_ranges.end(); allr_it++, nv_it++, et_it++) {
        // pack the entity type
      PACK_INT(buff_ptr, ((int)*et_it));

        // pack # ents
      PACK_INT(buff_ptr, (*allr_it).size());
      
        // pack the nodes per entity
      PACK_INT(buff_ptr, *nv_it);
      
        // pack the connectivity
      const MBEntityHandle *connect;
      int num_connect;
      MBEntityHandle *start_vec = (MBEntityHandle*)buff_ptr;
      for (MBRange::const_iterator rit = allr_it->begin(); rit != allr_it->end(); rit++) {
        result = mbImpl->get_connectivity(*rit, connect, num_connect);
        RRA("Failed to get connectivity.");
        assert(num_connect == *nv_it);
        PACK_EH(buff_ptr, &connect[0], num_connect);
      }

        // substitute destination handles
      result = get_remote_handles(store_remote_handles, start_vec, start_vec,
                                  *nv_it*(*allr_it).size(), to_proc,
                                  whole_range);
      RRA("Trouble getting remote handles when packing entities.");

        // pack the handles
      if (store_remote_handles)
        PACK_RANGE(buff_ptr, (*allr_it));

#ifdef DEBUG_PACKING
      std::cerr << "Packed " << (*allr_it).size() << " ents of type " 
                << MBCN::EntityTypeName(TYPE_FROM_HANDLE(*(allr_it->begin()))) << std::endl;
#endif      
      
    }

      // pack MBMAXTYPE to indicate end of ranges
    PACK_INT(buff_ptr, ((int)MBMAXTYPE));

    count += buff_ptr - orig_buff_ptr;
  }

  if (debug_packing) std::cerr << std::endl << "Done packing entities." << std::endl;

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_remote_handles(const bool store_remote_handles,
                                               MBEntityHandle *from_vec, 
                                               MBEntityHandle *to_vec_tmp,
                                               int num_ents, int to_proc,
                                               const MBRange &new_ents) 
{
    // NOTE: THIS IMPLEMENTATION IS JUST LIKE THE RANGE-BASED VERSION, NO REUSE
    // AT THIS TIME, SO IF YOU FIX A BUG IN THIS VERSION, IT MAY BE IN THE
    // OTHER VERSION TOO!!!
  if (0 == num_ents) return MB_SUCCESS;
  
    // use a local destination ptr in case we're doing an in-place copy
  std::vector<MBEntityHandle> tmp_vector;
  MBEntityHandle *to_vec = to_vec_tmp;
  if (to_vec == from_vec) {
    tmp_vector.resize(num_ents);
    to_vec = &tmp_vector[0];
  }

  if (!store_remote_handles) {
    int err;
      // in this case, substitute position in new_ents list
    for (int i = 0; i < num_ents; i++) {
      int ind = new_ents.index(from_vec[i]);
      to_vec[i] = CREATE_HANDLE(MBMAXTYPE, ind, err);
      assert(to_vec[i] != 0 && !err && -1 != ind);
    }
  }
  else {
    MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
    MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                              sharedh_tag, sharedhs_tag, pstatus_tag);
  
      // get single-proc destination handles and shared procs
    std::vector<int> sharing_procs(num_ents);
    result = mbImpl->tag_get_data(sharedh_tag, from_vec, num_ents,
                                  to_vec);
    RRA("Failed to get shared handle tag for remote_handles.");
    result = mbImpl->tag_get_data(sharedp_tag, from_vec, num_ents, &sharing_procs[0]);
    RRA("Failed to get sharing proc tag in remote_handles.");
    for (int j = 0; j < num_ents; j++) {
      if (to_vec[j] && sharing_procs[j] != to_proc)
        to_vec[j] = 0;
    }
    
    MBEntityHandle tmp_handles[MAX_SHARING_PROCS];
    int tmp_procs[MAX_SHARING_PROCS];
    int i;
      // go through results, and for 0-valued ones, look for multiple shared proc
    MBEntityHandle *tmp_eh;
    for (tmp_eh = to_vec, i = 0; i < num_ents; i++) {
      if (!to_vec[i]) {
        result = mbImpl->tag_get_data(sharedps_tag, from_vec+i, 1, tmp_procs);
        if (MB_SUCCESS == result) {
          for (int j = 0; j < MAX_SHARING_PROCS; j++) {
            if (-1 == tmp_procs[j]) break;
            else if (tmp_procs[j] == to_proc) {
              result = mbImpl->tag_get_data(sharedhs_tag, from_vec+i, 1, tmp_handles);
              RRA("Trouble getting sharedhs tag.");
              to_vec[i] = tmp_handles[j];
              assert(to_vec[i]);
              break;
            }
          }
        }
        if (!to_vec[i]) {
          int j = new_ents.index(from_vec[i]);
          if (-1 == j) {
            result = MB_FAILURE;
            RRA("Failed to find new entity in send list.");
          }
          int err;
          to_vec[i] = CREATE_HANDLE(MBMAXTYPE, j, err);
          if (err) {
            result = MB_FAILURE;
            RRA("Failed to create handle in remote_handles.");
          }
        }
      }
    }
  }
  
    // memcpy over results if from_vec and to_vec are the same
  if (to_vec_tmp == from_vec) 
    memcpy(from_vec, to_vec, num_ents * sizeof(MBEntityHandle));
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_remote_handles(const bool store_remote_handles,
                                               const MBRange &from_range, 
                                               MBEntityHandle *to_vec,
                                               int to_proc,
                                               const MBRange &new_ents) 
{
    // NOTE: THIS IMPLEMENTATION IS JUST LIKE THE VECTOR-BASED VERSION, NO REUSE
    // AT THIS TIME, SO IF YOU FIX A BUG IN THIS VERSION, IT MAY BE IN THE
    // OTHER VERSION TOO!!!
  if (from_range.empty()) return MB_SUCCESS;
  
  if (!store_remote_handles) {
    int err;
      // in this case, substitute position in new_ents list
    MBRange::iterator rit;
    unsigned int i;
    for (rit = from_range.begin(), i = 0; rit != from_range.end(); rit++, i++) {
      int ind = new_ents.index(*rit);
      to_vec[i] = CREATE_HANDLE(MBMAXTYPE, ind, err);
      assert(to_vec[i] != 0 && !err && -1 != ind);
    }
  }
  else {
    MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
    MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                              sharedh_tag, sharedhs_tag, pstatus_tag);
  
      // get single-proc destination handles and shared procs
    std::vector<int> sharing_procs(from_range.size());
    result = mbImpl->tag_get_data(sharedh_tag, from_range, to_vec);
    RRA("Failed to get shared handle tag for remote_handles.");
    result = mbImpl->tag_get_data(sharedp_tag, from_range, &sharing_procs[0]);
    RRA("Failed to get sharing proc tag in remote_handles.");
    for (unsigned int j = 0; j < from_range.size(); j++) {
      if (to_vec[j] && sharing_procs[j] != to_proc)
        to_vec[j] = 0;
    }
    
    MBEntityHandle tmp_handles[MAX_SHARING_PROCS];
    int tmp_procs[MAX_SHARING_PROCS];
      // go through results, and for 0-valued ones, look for multiple shared proc
    MBRange::iterator rit;
    unsigned int i;
    for (rit = from_range.begin(), i = 0; rit != from_range.end(); rit++, i++) {
      if (!to_vec[i]) {
        result = mbImpl->tag_get_data(sharedhs_tag, &(*rit), 1, tmp_handles);
        if (MB_SUCCESS == result) {
          result = mbImpl->tag_get_data(sharedps_tag, &(*rit), 1, tmp_procs);
          RRA("Trouble getting sharedps tag.");
          for (int j = 0; j < MAX_SHARING_PROCS; j++)
            if (tmp_procs[j] == to_proc) {
              to_vec[i] = tmp_handles[j];
              break;
            }
        }
      
        if (!to_vec[i]) {
          int j = new_ents.index(*rit);
          if (-1 == j) {
            result = MB_FAILURE;
            RRA("Failed to find new entity in send list.");
          }
          int err;
          to_vec[i] = CREATE_HANDLE(MBMAXTYPE, j, err);
          if (err) {
            result = MB_FAILURE;
            RRA("Failed to create handle in remote_handles.");
          }
        }
      }
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_remote_handles(const bool store_remote_handles,
                                               const MBRange &from_range, 
                                               MBRange &to_range,
                                               int to_proc,
                                               const MBRange &new_ents) 
{
  std::vector<MBEntityHandle> to_vector(from_range.size());

  MBErrorCode result =
    get_remote_handles(store_remote_handles, from_range, &to_vector[0],
                       to_proc, new_ents);
  RRA("Trouble getting remote handles.");
  std::copy(to_vector.begin(), to_vector.end(), mb_range_inserter(to_range));
  return result;
}

MBErrorCode MBParallelComm::unpack_entities(unsigned char *&buff_ptr,
                                            const bool store_remote_handles,
                                            const int from_proc,
                                            MBRange &entities) 
{
  MBErrorCode result;
  bool done = false;
  MBReadUtilIface *ru = NULL;
  result = mbImpl->query_interface(std::string("MBReadUtilIface"), 
                                   reinterpret_cast<void**>(&ru));
  RRA("Failed to get MBReadUtilIface.");
  
  
  while (!done) {
    MBEntityType this_type;
    MBEntityHandle actual_start;
    UNPACK_INT(buff_ptr, this_type);
    assert(this_type >= MBVERTEX && 
           (this_type == MBMAXTYPE || this_type < MBENTITYSET));

      // MBMAXTYPE signifies end of entities data
    if (MBMAXTYPE == this_type) break;
    
      // get the number of ents
    int num_ents;
    UNPACK_INT(buff_ptr, num_ents);
    
    if (MBVERTEX == this_type) {
        // unpack coords
      if (num_ents) {
        std::vector<double*> coords(3);
        result = ru->get_node_arrays(3, num_ents, 0, proc_config().proc_rank(), 
                                     actual_start, coords);
        RRA("Failed to allocate node arrays.");

        entities.insert(actual_start, actual_start+num_ents-1);
      
          // unpack the buffer data directly into coords
        for (int i = 0; i < 3; i++) 
          memcpy(coords[i], buff_ptr+i*num_ents*sizeof(double), 
                 num_ents*sizeof(double));
        buff_ptr += 3*num_ents * sizeof(double);
        UPC(3*num_ents, " doubles");

          // unpack source handles
        if (store_remote_handles) {
          MBRange this_range;
          this_range.insert(actual_start, actual_start+num_ents-1);
          MBRange dum_range;
          UNPACK_RANGE(buff_ptr, dum_range);
          result = set_remote_data(this_range, dum_range, from_proc);
          RRA("Couldn't set sharing data");
        }
      }
    }

    else {
      int verts_per_entity;
      
        // unpack the nodes per entity
      UNPACK_INT(buff_ptr, verts_per_entity);
      
      MBEntityHandle *connect;
      result = ru->get_element_array(num_ents, verts_per_entity, this_type,
                                     0, proc_config().proc_rank(), actual_start,
                                     connect);
      RRA("Failed to allocate element arrays.");

        // unpack the connectivity
      UNPACK_EH(buff_ptr, connect, (num_ents*verts_per_entity));
      entities.insert(actual_start, actual_start+num_ents-1);

        // convert to local handles
      result = get_local_handles(connect, num_ents*verts_per_entity,
                                 entities);
      RRA("Couldn't get local handles.");

      if (store_remote_handles) {
          // unpack source handles
        MBRange this_range;
        this_range.insert(actual_start, actual_start+num_ents-1);
        MBRange dum_range;
        UNPACK_RANGE(buff_ptr, dum_range);
        result = set_remote_data(this_range, dum_range, from_proc);
        RRA("Couldn't set sharing data");
      }
    }

#ifdef DEBUG_PACKING
      std::cerr << "Unpacked " << num_ents << " ents of type " 
                << MBCN::EntityTypeName(TYPE_FROM_HANDLE(actual_start)) << std::endl;
#endif      

  }
  
  
  if (debug_packing) std::cerr << std::endl << "Done unpacking entities." << std::endl;

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::unpack_iface_entities(unsigned char *&buff_ptr,
                                                  const int from_proc, 
                                                  std::vector<MBEntityHandle> &recd_ents) 
{
    // for all the entities in the received buffer; for each, save
    // entities in this instance which match connectivity, or zero if none found
  MBErrorCode result;
  bool done = false;
  std::vector<MBEntityHandle> connect;
  
  while (!done) {
    MBEntityType this_type;
    UNPACK_INT(buff_ptr, this_type);
    assert(this_type > MBVERTEX && 
           (this_type == MBMAXTYPE || this_type < MBENTITYSET));

      // MBMAXTYPE signifies end of entities data
    if (MBMAXTYPE == this_type) break;
    
      // get the number of ents
    int num_ents;
    UNPACK_INT(buff_ptr, num_ents);
    
    int verts_per_entity;
      
      // unpack the verts per entity
    UNPACK_INT(buff_ptr, verts_per_entity);

    connect.resize(verts_per_entity * num_ents);
    
      // unpack the connectivity
    UNPACK_EH(buff_ptr, &connect[0], (num_ents*verts_per_entity));

      // unpack source handles
    MBRange source_range, tmp_range;
    UNPACK_RANGE(buff_ptr, source_range);

      // save matching local entities, or 0 if none found
    for (int i = 0; i < num_ents; i++) {
        // test for existing entity
      result = mbImpl->get_adjacencies(&connect[0]+i*verts_per_entity, verts_per_entity, 
                                       MBCN::Dimension(this_type), false, tmp_range);
      if (MB_MULTIPLE_ENTITIES_FOUND == result) {
        RRA("Multiple entities found.");
      }
      if (!tmp_range.empty()) {
          // found a corresponding entity - save target
        assert(1 == tmp_range.size());
        recd_ents.push_back(*tmp_range.begin());
        tmp_range.clear();
      }
      else
        recd_ents.push_back(0);
    }

#ifdef DEBUG_PACKING
    std::cerr << "Unpacked " << num_ents << " ents of type " 
              << MBCN::EntityTypeName(TYPE_FROM_HANDLE(actual_start)) << std::endl;
#endif      

  }
  
  
  if (debug_packing) std::cerr << std::endl << "Done unpacking iface entities." << std::endl;

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_local_handles(const MBRange &remote_handles,
                                              MBRange &local_handles,
                                              const MBRange &new_ents) 
{
  std::vector<MBEntityHandle> rh_vec;
  rh_vec.reserve(remote_handles.size());
  std::copy(remote_handles.begin(), remote_handles.end(), std::back_inserter(rh_vec));
  MBErrorCode result = get_local_handles(&rh_vec[0], remote_handles.size(), new_ents);
  std::copy(rh_vec.begin(), rh_vec.end(), mb_range_inserter(local_handles));
  return result;
}
  
MBErrorCode MBParallelComm::get_local_handles(MBEntityHandle *from_vec, 
                                              int num_ents,
                                              const MBRange &new_ents) 
{
  for (int i = 0; i < num_ents; i++) {
    if (TYPE_FROM_HANDLE(from_vec[i]) == MBMAXTYPE) {
      assert(ID_FROM_HANDLE(from_vec[i]) < (int) new_ents.size());
      from_vec[i] = new_ents[ID_FROM_HANDLE(from_vec[i])];
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::set_remote_data(MBRange &local_range,
                                            MBRange &remote_range,
                                            int other_proc) 
{
    // NOTE: THIS IMPLEMENTATION IS JUST LIKE THE VECTOR-BASED VERSION, NO REUSE
    // AT THIS TIME, SO IF YOU FIX A BUG IN THIS VERSION, IT MAY BE IN THE
    // OTHER VERSION TOO!!!

  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag, pstatus_tag);

    // get remote procs tag, if any
  MBRange tmp_range, tmp_local_range = local_range;
  std::vector<int> remote_proc(local_range.size());
  std::vector<int> remote_procs(MAX_SHARING_PROCS);
  std::vector<MBEntityHandle> remote_handle(local_range.size());
  std::vector<MBEntityHandle> remote_handles(MAX_SHARING_PROCS);
  std::fill(remote_procs.begin(), remote_procs.end(), -1);
  std::fill(remote_handles.begin(), remote_handles.end(), 0);
  result = mbImpl->tag_get_data(sharedp_tag, local_range,
                                &remote_proc[0]);
  RRA("Couldn't get sharedp tag (range).");
  result = mbImpl->tag_get_data(sharedh_tag, local_range,
                                &remote_handle[0]);
  RRA("Couldn't get sharedh tag (range).");
  MBRange::iterator rit, rit2;
  int i = 0;

    // for each pair of local/remote handles:
  for (rit = tmp_local_range.begin(), rit2 = remote_range.begin(); 
       rit != tmp_local_range.end(); rit++, rit2++, i++) {

      // get existing remote proc(s), handle(s) for this local handle
    remote_procs[0] = remote_proc[i];
    if (-1 != remote_procs[0]) remote_handles[0] = remote_handle[i];
    else {
      result = mbImpl->tag_get_data(sharedps_tag, &(*rit), 1,
                                    &remote_procs[0]);
      if (MB_SUCCESS == result) {
        result = mbImpl->tag_get_data(sharedhs_tag, &(*rit), 1,
                                      &remote_handles[0]);
        RRA("Couldn't get sharedhs tag (range).");
      }
    }

      // now insert other_proc, handle into these, and set the tags
    if (-1 == remote_procs[0]) {
      remote_procs[0] = other_proc;
      remote_handles[0] = *rit2;
      result = mbImpl->tag_set_data(sharedp_tag, &(*rit), 1, &remote_procs[0]);
      RRA("Couldn't set sharedp tag (range)");
      result = mbImpl->tag_set_data(sharedh_tag, &(*rit), 1, &remote_handles[0]);
      RRA("Couldn't set sharedh tag (range)");
      remote_procs[0] = -1;
      remote_handles[0] = 0;
    }
    else {
      std::vector<int>::iterator vit = remote_procs.begin();
      std::vector<MBEntityHandle>::iterator vit2 = remote_handles.begin();
      while (-1 != *vit && vit != remote_procs.end() && other_proc > *vit) {
        vit++; vit2++;
      }
      assert(*vit != other_proc);
      remote_procs.insert(vit, other_proc);
      remote_handles.insert(vit2, *rit2);
        // need to pop 'cuz we're using an insert above
      remote_procs.pop_back();
      remote_handles.pop_back();
      result = mbImpl->tag_set_data(sharedps_tag, &(*rit), 1, &remote_procs[0]);
      RRA("Couldn't set sharedps tag (range)");
      result = mbImpl->tag_set_data(sharedhs_tag, &(*rit), 1, &remote_handles[0]);
      RRA("Couldn't set sharedhs tag (range)");
      std::fill(remote_procs.begin(), remote_procs.end(), -1);
      std::fill(remote_handles.begin(), remote_handles.end(), 0);
    }
  }

    // also update shared flag for these ents
  unsigned int *shared_flags = (unsigned int*) &remote_proc[0];
  result = mbImpl->tag_get_data(pstatus_tag, local_range, shared_flags);
  RRA("Couldn't get pstatus tag (range)");
  for (unsigned int i = 0; i < local_range.size(); i++)
    shared_flags[i] |= PSTATUS_SHARED;
  result = mbImpl->tag_set_data(pstatus_tag, local_range, shared_flags);
  RRA("Couldn't set pstatus tag (range)");
  
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::set_remote_data(MBEntityHandle *local_ents,
                                            MBEntityHandle *remote_ents,
                                            int num_ents,
                                            int other_proc) 
{
    // NOTE: THIS IMPLEMENTATION IS JUST LIKE THE RANGE-BASED VERSION, NO REUSE
    // AT THIS TIME, SO IF YOU FIX A BUG IN THIS VERSION, IT MAY BE IN THE
    // OTHER VERSION TOO!!!
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag, pstatus_tag);

    // get remote procs tag, if any
  MBRange tmp_range;
  std::vector<int> remote_proc(num_ents);
  int remote_procs[MAX_SHARING_PROCS];
  std::vector<MBEntityHandle> remote_handle(num_ents);
  MBEntityHandle remote_handles[MAX_SHARING_PROCS];
  std::vector<unsigned char> pstatus_vals(num_ents);
  result = mbImpl->tag_get_data(sharedp_tag, local_ents, num_ents,
                                &remote_proc[0]);
  RRA("Couldn't get sharedp tag (vector)");
  result = mbImpl->tag_get_data(sharedh_tag, local_ents, num_ents,
                                &remote_handle[0]);
  RRA("Couldn't get sharedh tag (vector)");
  result = mbImpl->tag_get_data(pstatus_tag, local_ents, num_ents,
                                &pstatus_vals[0]);
  RRA("Couldn't get sharedh tag (vector)");

    // for each local/remote handle pair
  for (int i = 0; i != num_ents; i++) {

      // get existing remote proc(s), handle(s) for this local handle
    if (-1 != remote_proc[i]) {
      remote_procs[0] = remote_proc[i];
      remote_handles[0] = remote_handle[i];
      std::fill(remote_procs+1, remote_procs+MAX_SHARING_PROCS, -1);
      std::fill(remote_handles+1, remote_handles+MAX_SHARING_PROCS, 0);
    }
    else {
      result = mbImpl->tag_get_data(sharedps_tag, local_ents+i, 1,
                                    remote_procs);
      if (MB_SUCCESS == result) {
        result = mbImpl->tag_get_data(sharedhs_tag, local_ents+i, 1,
                                      remote_handles);
        RRA("Couldn't get sharedhs tag (vector)");
      }
    }

      // now either insert other_proc, handle into these, or remove if
      // remote handle is 0
    if (0 == remote_ents[i]) {
      result = rmv_remote_proc(local_ents[i], remote_procs, remote_handles,
                               other_proc);
      RRA(" ");
    }
    else {
      result = add_remote_proc(local_ents[i], remote_procs, remote_handles,
                               other_proc, remote_ents[i]);
      RRA(" ");
    }
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::rmv_remote_proc(MBEntityHandle ent,
                                            int *remote_procs,
                                            MBEntityHandle *remote_hs,
                                            int remote_proc) 
{
  int i = std::find(remote_procs, remote_procs+MAX_SHARING_PROCS, -1) - remote_procs;
  int j = std::find(remote_procs, remote_procs+i, remote_proc) - remote_procs;
  MBErrorCode result;

    // remote_proc must be in this list
  assert(j != i);
  
  for (; j != i; j++) {
    remote_procs[j] = remote_procs[j+1];
    remote_hs[j] = remote_hs[j+1];
  }

  if (i == 1) {
      // only had one proc, and now have none; have to set the values,
      // since they're sparse tags
    result = mbImpl->tag_set_data(sharedp_tag(), &ent, 1, remote_procs);
    result = mbImpl->tag_set_data(sharedh_tag(), &ent, 1, remote_hs);
      // if it's not shared, it's also not a ghost, interface, or !owned entity
    unsigned char not_shared = 0;
    result = mbImpl->tag_set_data(pstatus_tag(), &ent, 1, &not_shared);
  }
  else if (i == 2) {
      // went from 2 to 1, need to unset 1 and set the other
    result = mbImpl->tag_set_data(sharedp_tag(), &ent, 1, remote_procs);
    RRA("Couldn't set sharedp tag");
    result = mbImpl->tag_set_data(sharedh_tag(), &ent, 1, remote_hs);
    RRA("Couldn't set sharedh tag");
    result = mbImpl->tag_delete_data(sharedps_tag(), &ent, 1);
    RRA("Couldn't remove sharedps tag");
    result = mbImpl->tag_delete_data(sharedhs_tag(), &ent, 1);
    RRA("Couldn't remove sharedhs tag");
  }
  else {
    result = mbImpl->tag_set_data(sharedps_tag(), &ent, 1, remote_procs);
    RRA("Couldn't set sharedps tag");
    result = mbImpl->tag_set_data(sharedhs_tag(), &ent, 1, remote_hs);
    RRA("Couldn't set sharedhs tag");
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::add_remote_proc(MBEntityHandle ent,
                                            int *remote_procs,
                                            MBEntityHandle *remote_hs,
                                            int remote_proc,
                                            MBEntityHandle remote_handle) 
{
  int i = 0;
  bool more = (remote_procs[0] != -1);
  
  while (-1 != remote_procs[i] && i < MAX_SHARING_PROCS && 
         remote_proc > remote_procs[i]) i++;

  if (remote_procs[i] == remote_proc) return MB_SUCCESS;
  
  assert(i < MAX_SHARING_PROCS-1);
  if (-1 != remote_procs[i]) {
    for (int j = MAX_SHARING_PROCS-1; j > i; j--) {
      remote_procs[j] = remote_procs[j-1];
      remote_hs[j] = remote_hs[j-1];
    }
  }
  remote_procs[i] = remote_proc;
  remote_hs[i] = remote_handle;

  MBErrorCode result = mbImpl->tag_set_data((more ? sharedps_tag() : sharedp_tag()), 
                                            &ent, 1, remote_procs);
  RRA("Couldn't set sharedps tag");
  result = mbImpl->tag_set_data((more ? sharedhs_tag() : sharedh_tag()), &ent, 1, remote_hs);
  RRA("Couldn't set sharedhs tag");

    // if we went from 1 to 2, need to unset sharedp_tag
  if (1 == i) {
    remote_proc = -1;
    result = mbImpl->tag_set_data(sharedp_tag(), &ent, 1, &remote_proc);
    RRA("Couldn't set sharedp tag");
    remote_handle = 0;
    result = mbImpl->tag_set_data(sharedh_tag(), &ent, 1, &remote_handle);
    RRA("Couldn't set sharedhs tag");
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_range_map(MBRange &key_range, MBEntityHandle val_start,
                                           HandleMap &handle_map) 
{
  for (MBRange::const_pair_iterator key_it = key_range.const_pair_begin(); 
       key_it != key_range.const_pair_end(); key_it++) {
    int tmp_num = (*key_it).second - (*key_it).first + 1;
    handle_map.insert((*key_it).first, val_start, tmp_num);
    val_start += tmp_num;
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_sets(MBRange &entities,
                                      MBRange::const_iterator &start_rit,
                                      MBRange &whole_range,
                                      unsigned char *&buff_ptr,
                                      int &count,
                                      const bool just_count,
                                      const bool store_remote_handles,
                                      const int to_proc,
                                      MBRange &set_range,
                                      std::vector<MBRange> &set_ranges,
                                      std::vector<int> &set_sizes,
                                      std::vector<unsigned int> &options_vec)
{
    // SETS:
    // . #sets
    // . for each set:
    //   - options[#sets] (unsigned int)
    //   - if (unordered) set range 
    //   - else if ordered
    //     . #ents in set
    //     . handles[#ents]
    //   - #parents
    //   - if (#parents) handles[#parents]
    //   - #children
    //   - if (#children) handles[#children]
  
    // now the sets; assume any sets the application wants to pass are in the entities list
  unsigned char *orig_buff_ptr = buff_ptr;
  MBErrorCode result;

  if (just_count) {
    int ranges_size = 0, vecs_size = 0, tot_parch = 0;
    
    for (; start_rit != entities.end(); start_rit++) {
      set_range.insert(*start_rit);
    
      unsigned int options;
      result = mbImpl->get_meshset_options(*start_rit, options);
      RRA("Failed to get meshset options.");
      options_vec.push_back(options);

      if (options & MESHSET_SET) {
          // range-based set; count the subranges
        set_ranges.push_back(MBRange());
        result = mbImpl->get_entities_by_handle(*start_rit, *set_ranges.rbegin());
        RRA("Failed to get set entities.");

          // set range
        ranges_size += 2 * sizeof(MBEntityHandle) * num_subranges(*set_ranges.rbegin()) + 
          sizeof(int);
      }
      else if (options & MESHSET_ORDERED) {
          // just get the number of entities in the set
        int num_ents;
        result = mbImpl->get_number_entities_by_handle(*start_rit, num_ents);
        RRA("Failed to get number entities in ordered set.");
        set_sizes.push_back(num_ents);

          // set vec
        vecs_size += sizeof(MBEntityHandle) * num_ents + sizeof(int);
      }

        // get numbers of parents/children
      int num_par, num_ch;
      result = mbImpl->num_child_meshsets(*start_rit, &num_ch);
      RRA("Failed to get num children.");

      result = mbImpl->num_parent_meshsets(*start_rit, &num_par);
      RRA("Failed to get num parents.");

      tot_parch += num_ch + num_par;
      
    }

      // num of sets
    count += sizeof(int);
    
      // options
    count += set_range.size() * sizeof(unsigned int);

      // range, vector sizes
    count += ranges_size + vecs_size;
    
      // set children, parents
    count += 2 * set_range.size() * sizeof(int) + tot_parch * sizeof(MBEntityHandle);

        // set handles
    if (!set_range.empty() && store_remote_handles)
      count += sizeof(int) + 2*sizeof(MBEntityHandle)*num_subranges(set_range);

    whole_range.merge(set_range);
  }
  else {
    
      // set handle range
    PACK_INT(buff_ptr, set_range.size());

      // option values
    if (!set_range.empty())
      PACK_VOID(buff_ptr, &options_vec[0], set_range.size()*sizeof(unsigned int));

    std::vector<unsigned int>::const_iterator opt_it = options_vec.begin();
    std::vector<MBRange>::const_iterator rit = set_ranges.begin();
    std::vector<int>::const_iterator mem_it = set_sizes.begin();
    static std::vector<MBEntityHandle> members;

    for (MBRange::const_iterator set_it = set_range.begin(); set_it != set_range.end(); 
         set_it++, opt_it++) {
      if ((*opt_it) & MESHSET_SET) {
          // pack entities as a range
        MBRange dum_range;
        result = get_remote_handles(store_remote_handles,
                                    (*rit), dum_range, to_proc,
                                    whole_range);
        RRA("Trouble getting remote handles for unordered set contents.");
        PACK_RANGE(buff_ptr, dum_range);
        rit++;
      }
      else if ((*opt_it) & MESHSET_ORDERED) {
          // pack entities as vector, with length
        PACK_INT(buff_ptr, *mem_it);
        members.clear();
        result = mbImpl->get_entities_by_handle(*set_it, members);
        RRA("Failed to get set entities.");
        result = get_remote_handles(store_remote_handles,
                                    &members[0], &members[0], 
                                    members.size(), to_proc,
                                    whole_range);
        RRA("Trouble getting remote handles for ordered set contents.");
        PACK_EH(buff_ptr, &members[0], *mem_it);
        mem_it++;
      }
    }
    
      // pack numbers of parents/children
    unsigned int tot_pch = 0;
    for (MBRange::const_iterator set_it = set_range.begin(); set_it != set_range.end(); 
         set_it++) {
        // pack parents
      int num_pch;
      result = mbImpl->num_parent_meshsets(*set_it, &num_pch);
      RRA("Failed to get num parents.");
      PACK_INT(buff_ptr, num_pch);
      tot_pch += num_pch;
      result = mbImpl->num_child_meshsets(*set_it, &num_pch);
      RRA("Failed to get num children.");
      PACK_INT(buff_ptr, num_pch);
      tot_pch += num_pch;
    }

      // now pack actual parents/children
    members.clear();
    members.reserve(tot_pch);
    std::vector<MBEntityHandle> tmp_pch;
    for (MBRange::const_iterator set_it = set_range.begin(); set_it != set_range.end(); 
         set_it++) {
      result = mbImpl->get_parent_meshsets(*set_it, tmp_pch);
      RRA("Failed to get parents.");
      std::copy(tmp_pch.begin(), tmp_pch.end(), std::back_inserter(members));
      tmp_pch.clear();
      result = mbImpl->get_child_meshsets(*set_it, tmp_pch);
      RRA("Failed to get children.");
      std::copy(tmp_pch.begin(), tmp_pch.end(), std::back_inserter(members));
      tmp_pch.clear();
    }
    assert(members.size() == tot_pch);
    if (!members.empty()) {
      result = get_remote_handles(store_remote_handles,
                                  &members[0], &members[0], 
                                  members.size(), to_proc,
                                  whole_range);
      RRA("Trouble getting remote handles for set parent/child sets.");
#ifndef NDEBUG
        // check that all handles are either sets or maxtype
      for (unsigned int __j = 0; __j < members.size(); __j++)
        assert((TYPE_FROM_HANDLE(members[__j]) == MBMAXTYPE &&
                ID_FROM_HANDLE(members[__j]) < (int)whole_range.size()) ||
               TYPE_FROM_HANDLE(members[__j]) == MBENTITYSET);
#endif        
      PACK_EH(buff_ptr, &members[0], members.size());
    }
    
    // pack the handles
    if (store_remote_handles && !set_range.empty())
      PACK_RANGE(buff_ptr, set_range);

    count += buff_ptr - orig_buff_ptr;
  }

  if (debug_packing) std::cerr << std::endl << "Done packing sets." << std::endl;

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::unpack_sets(unsigned char *&buff_ptr,
                                        MBRange &entities,
                                        const bool store_remote_handles,
                                        const int from_proc)
{
  
    // now the sets; assume any sets the application wants to pass are in the entities list
  MBErrorCode result;

  MBRange new_sets;
  int num_sets;
  UNPACK_INT(buff_ptr, num_sets);

  if (!num_sets) return MB_SUCCESS;
         
  std::vector<MBEntityHandle> members;
  int num_ents;
  std::vector<unsigned int> options_vec(num_sets);
      // option value
  if (num_sets)
    UNPACK_VOID(buff_ptr, &options_vec[0], num_sets*sizeof(unsigned int));

    // create sets
  int i;
  MBRange::const_iterator rit;
  for (i = 0; i < num_sets; i++) {
    
      // create the set
    MBEntityHandle set_handle;
    result = mbImpl->create_meshset(options_vec[i], set_handle);
    RRA("Failed to create set in unpack.");
    new_sets.insert(set_handle);
  }

    // make sure new sets handles are contiguous
  assert(*new_sets.rbegin() - *new_sets.begin() + 1 == new_sets.size());

  entities.merge(new_sets);
  
  for (rit = new_sets.begin(), i = 0; rit != new_sets.end(); rit++, i++) {
    if (options_vec[i] & MESHSET_SET) {
        // unpack entities as a range
      MBRange set_range, tmp_range;
      UNPACK_RANGE(buff_ptr, tmp_range);
      result = get_local_handles(tmp_range, set_range, entities);      
      RRA("Failed to get local handles for unordered set contents.");
      result = mbImpl->add_entities(*rit, set_range);
      RRA("Failed to add ents to unordered set in unpack.");
    }
    else if (options_vec[i] & MESHSET_ORDERED) {
        // unpack entities as vector, with length
      UNPACK_INT(buff_ptr, num_ents);
      members.resize(num_ents);
      if (num_ents) UNPACK_EH(buff_ptr, &members[0], num_ents);
      result = get_local_handles(&members[0], num_ents, entities);
      RRA("Failed to get local handles for ordered set contents.");
      result = mbImpl->add_entities(*rit, &members[0], num_ents);
      RRA("Failed to add ents to ordered set in unpack.");
    }
  }

  std::vector<int> num_pch(2*new_sets.size());
  std::vector<int>::iterator vit;
  int tot_pch = 0;
  for (vit = num_pch.begin(); vit != num_pch.end(); vit++) {
    UNPACK_INT(buff_ptr, *vit);
    tot_pch += *vit;
  }
  
  members.resize(tot_pch);
  UNPACK_EH(buff_ptr, &members[0], tot_pch);
  result = get_local_handles(&members[0], tot_pch, entities);
  RRA("Couldn't get local handle for parent/child sets.");

  int num = 0;
  MBEntityHandle *mem_ptr = &members[0];
  for (rit = new_sets.begin(); rit != new_sets.end(); rit++) {
      // unpack parents/children
    int num_par = num_pch[num++], num_child = num_pch[num++];
    if (num_par+num_child) {
      for (i = 0; i < num_par; i++) {
        assert(0 != mem_ptr[i]);
        result = mbImpl->add_parent_meshset(*rit, mem_ptr[i]);
        RRA("Failed to add parent to set in unpack.");
      }
      mem_ptr += num_par;
      for (i = 0; i < num_child; i++) {
        assert(0 != mem_ptr[i]);
        result = mbImpl->add_child_meshset(*rit, mem_ptr[i]);
        RRA("Failed to add child to set in unpack.");
      }
      mem_ptr += num_child;
    }
  }

    // unpack source handles
  MBRange dum_range;
  if (store_remote_handles && !new_sets.empty()) {
    UNPACK_RANGE(buff_ptr, dum_range);
    result = set_remote_data(new_sets, dum_range, from_proc);
    RRA("Couldn't set sharing data for sets");
  }

  if (debug_packing) std::cerr << std::endl << "Done unpacking sets." << std::endl;

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_adjacencies(MBRange &entities,
                                             MBRange::const_iterator &start_rit,
                                             MBRange &whole_range,
                                             unsigned char *&buff_ptr,
                                             int &count,
                                             const bool just_count,
                                             const bool store_handles,
                                             const int to_proc)
{
  return MB_FAILURE;
}

MBErrorCode MBParallelComm::unpack_adjacencies(unsigned char *&buff_ptr,
                                               MBRange &entities,
                                               const bool store_handles,
                                               const int from_proc)
{
  return MB_FAILURE;
}

MBErrorCode MBParallelComm::pack_tags(MBRange &entities,
                                      MBRange::const_iterator &start_rit,
                                      MBRange &whole_range,
                                      unsigned char *&buff_ptr,
                                      int &count,
                                      const bool just_count,
                                      const bool store_remote_handles,
                                      const int to_proc,
                                      std::vector<MBTag> &all_tags,
                                      std::vector<MBRange> &tag_ranges)
{
    // tags
    // get all the tags
    // for dense tags, compute size assuming all entities have that tag
    // for sparse tags, get number of entities w/ that tag to compute size

  unsigned char *orig_buff_ptr = buff_ptr;
  MBErrorCode result;
  std::vector<int> var_len_sizes;
  std::vector<const void*> var_len_values;

  if (just_count) {

    std::vector<MBTag> tmp_tags;
    result = tagServer->get_tags(tmp_tags);
    RRA("Failed to get tags in pack_tags.");

    for (std::vector<MBTag>::iterator tag_it = tmp_tags.begin(); tag_it != tmp_tags.end(); tag_it++) {
      const TagInfo *tinfo = tagServer->get_tag_info(*tag_it);
      MBRange tmp_range;
      result = tagServer->get_entities(*tag_it, tmp_range);
      RRA("Failed to get entities for tag in pack_tags.");
      tmp_range = tmp_range.intersect(whole_range);
      
      if (tmp_range.empty())
        continue;

        // ok, we'll be sending this tag
      all_tags.push_back(*tag_it);

        // default value
      count += sizeof(int);
      if (NULL != tinfo->default_value()) 
        count += tinfo->default_value_size();
      
        // size, type, data type
      count += 3*sizeof(int);
      
        // name
      count += sizeof(int);
      count += tinfo->get_name().size();

      if (!tmp_range.empty()) {
        tag_ranges.push_back(tmp_range);

        // range of tag
        count += sizeof(int) + tmp_range.size() * sizeof(MBEntityHandle);
      
        if (tinfo->get_size() == MB_VARIABLE_LENGTH) {
          const int num_ent = tmp_range.size();
            // send a tag size for each entity
          count += num_ent * sizeof(int);
            // send tag data for each entity
          var_len_sizes.resize( num_ent );
          var_len_values.resize( num_ent );
          result = tagServer->get_data( *tag_it, tmp_range, &var_len_values[0], 
                                        &var_len_sizes[0] );
          RRA("Failed to get lenghts of variable-length tag values.");
          count += std::accumulate( var_len_sizes.begin(), var_len_sizes.end(), 0 );
        }
        else {
            // tag data values for range or vector
          count += tmp_range.size() * tinfo->get_size();
        }
      }
    }
    
      // number of tags
    count += sizeof(int);
  }

  else {
    std::vector<unsigned char> tag_data;
    std::vector<MBRange>::const_iterator tr_it = tag_ranges.begin();

    PACK_INT(buff_ptr, all_tags.size());
    
    for (std::vector<MBTag>::const_iterator tag_it = all_tags.begin(); tag_it != all_tags.end(); tag_it++) {

      const TagInfo *tinfo = tagServer->get_tag_info(*tag_it);

        // size, type, data type
      PACK_INT(buff_ptr, tinfo->get_size());
      MBTagType this_type;
      result = mbImpl->tag_get_type(*tag_it, this_type);
      PACK_INT(buff_ptr, this_type);
      PACK_INT(buff_ptr, tinfo->get_data_type());
      
        // default value
      if (NULL == tinfo->default_value()) {
        PACK_INT(buff_ptr, 0);
      }
      else {
        PACK_INT(buff_ptr, tinfo->default_value_size());
        PACK_VOID(buff_ptr, tinfo->default_value(), tinfo->default_value_size());
      }
      
        // name
      PACK_INT(buff_ptr, tinfo->get_name().size() );
      PACK_VOID(buff_ptr, tinfo->get_name().c_str(), tinfo->get_name().size());
      
#ifdef DEBUG_PACKING
    std::cerr << "Packing tag " << tinfo->get_name() << std::endl;
#endif    
        // pack entities
      PACK_INT(buff_ptr, (*tr_it).size());
      result = get_remote_handles(store_remote_handles,
                                  (*tr_it), (MBEntityHandle*)buff_ptr, to_proc,
                                  whole_range);
#ifdef DEBUG_PACKING
      if (MB_SUCCESS != result) {
        std::cerr << "Trouble getting remote handles for tagged entities:" << std::endl;
        (*tr_it).print("  ");
      }
#else
      RRA("Trouble getting remote handles for tagged entities.");
#endif

      buff_ptr += (*tr_it).size() * sizeof(MBEntityHandle);

      const size_t num_ent = tr_it->size();
      if (tinfo->get_size() == MB_VARIABLE_LENGTH) {
        var_len_sizes.resize( num_ent, 0 );
        var_len_values.resize( num_ent, 0 );
        result = mbImpl->tag_get_data(*tag_it, *tr_it, &var_len_values[0], 
                                      &var_len_sizes[0] );
        RRA("Failed to get variable-length tag data in pack_tags.");
        PACK_INTS(buff_ptr, &var_len_sizes[0], num_ent);
        for (unsigned int i = 0; i < num_ent; ++i)
          PACK_VOID(buff_ptr, var_len_values[i], var_len_sizes[i]);
      }
      else {
        result = mbImpl->tag_get_data(*tag_it, *tr_it, buff_ptr);
        RRA("Failed to get tag data in pack_tags.");
        buff_ptr += num_ent * tinfo->get_size();
        PC(num_ent*tinfo->get_size(), " void");
      }
      tr_it++;
    }

    count += buff_ptr - orig_buff_ptr;
  }
  
  if (debug_packing) std::cerr << std::endl << "Done packing tags." << std::endl;

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::unpack_tags(unsigned char *&buff_ptr,
                                        MBRange &entities,
                                        const bool store_remote_handles,
                                        const int from_proc)
{
    // tags
    // get all the tags
    // for dense tags, compute size assuming all entities have that tag
    // for sparse tags, get number of entities w/ that tag to compute size

  MBErrorCode result;
  
  int num_tags;
  UNPACK_INT(buff_ptr, num_tags);
  std::vector<MBEntityHandle> tag_ents;
  std::vector<const void*> var_len_vals;
  std::vector<int> var_lengths;

  for (int i = 0; i < num_tags; i++) {
    
        // tag handle
    MBTag tag_handle;

      // size, data type
    int tag_size, tag_data_type, tag_type;
    UNPACK_INT(buff_ptr, tag_size);
    UNPACK_INT(buff_ptr, tag_type);
    UNPACK_INT(buff_ptr, tag_data_type);
      
      // default value
    int def_val_size;
    UNPACK_INT(buff_ptr, def_val_size);
    void *def_val_ptr = NULL;
    if (def_val_size) {
      def_val_ptr = buff_ptr;
      buff_ptr += def_val_size;
      UPC(tag_size, " void");
    }
    
      // name
    int name_len;
    UNPACK_INT(buff_ptr, name_len);
    std::string tag_name( reinterpret_cast<char*>(buff_ptr), name_len );
    buff_ptr += name_len;
    UPC(64, " chars");
#ifdef DEBUG_PACKING
    std::cerr << "Unpacking tag " << tag_name << std::endl;
#endif    

      // create the tag
    if (tag_size == MB_VARIABLE_LENGTH) 
      result = mbImpl->tag_create_variable_length( tag_name.c_str(), (MBTagType)tag_type,
                                                   (MBDataType)tag_data_type, tag_handle,
                                                   def_val_ptr, def_val_size );
    else
      result = mbImpl->tag_create(tag_name.c_str(), tag_size, (MBTagType) tag_type, 
                                  (MBDataType) tag_data_type, tag_handle,
                                  def_val_ptr);
    if (MB_ALREADY_ALLOCATED == result) {
        // already allocated tag, check to make sure it's the same size, type, etc.
      const TagInfo *tag_info = tagServer->get_tag_info(tag_name.c_str());
      MBTagType this_type;
      result = mbImpl->tag_get_type(tag_handle, this_type);
      if (tag_size != tag_info->get_size() ||
          tag_type != this_type ||
          tag_data_type != tag_info->get_data_type() ||
          (def_val_ptr && !tag_info->default_value()) ||
          (!def_val_ptr && tag_info->default_value())) {
        RRA("Didn't get correct tag info when unpacking tag.");
      }
    }
    else if (MB_SUCCESS != result) return result;

      // go through handle vec (in buffer) and convert to local handles in-place
    int num_ents;
    UNPACK_INT(buff_ptr, num_ents);
    MBEntityHandle *handle_vec = (MBEntityHandle*)buff_ptr;
    result = get_local_handles(handle_vec, num_ents, entities);
    RRA("Failed to get local handles for tagged entities.");
    buff_ptr += num_ents * sizeof(MBEntityHandle);

      // if it's a handle type, also convert tag vals in-place in buffer
    if (MB_TYPE_HANDLE == tag_type) {
      MBEntityHandle *val_vec = (MBEntityHandle*)buff_ptr;
      result = get_local_handles(val_vec, num_ents, entities);
      RRA("Failed to get local handles for tag vals.");
    }

    if (tag_size == MB_VARIABLE_LENGTH) {
        // Be careful of alignment here.  If the integers are aligned
        // in the buffer, we can use them directly.  Otherwise we must
        // copy them.
      const int* size_arr;
      if (((size_t)buff_ptr)%4) {
        var_lengths.resize( num_ents );
        memcpy( &var_lengths[0], buff_ptr, num_ents*sizeof(int) );
        size_arr = &var_lengths[0];
      }
      else {
        size_arr = reinterpret_cast<const int*>(buff_ptr);
      }
      buff_ptr += sizeof(int) * num_ents;
      UPC(sizeof(int) * num_ents, " void");
      
        // get pointers into buffer for each tag value
      var_len_vals.resize(num_ents);
      for (std::vector<MBEntityHandle>::size_type i = 0; 
           i < (std::vector<MBEntityHandle>::size_type) num_ents; ++i) {
        var_len_vals[i] = buff_ptr;
        buff_ptr += size_arr[i];
        UPC(size_arr[i], " void");
      }
      result = mbImpl->tag_set_data( tag_handle, handle_vec, num_ents,
                                     &var_len_vals[0], size_arr );
      RRA("Trouble setting tag data when unpacking variable-length tag.");
    }
    else {
      result = mbImpl->tag_set_data(tag_handle, handle_vec,
                                    num_ents, buff_ptr);
      RRA("Trouble setting range-based tag data when unpacking tag.");
      buff_ptr += num_ents * tag_size;
      UPC(num_ents * tag_size, " void");
    }
  }
  
  if (debug_packing) std::cerr << std::endl << "Done unpacking tags." << std::endl;

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::resolve_shared_ents(int resolve_dim,
                                                int shared_dim) 
{
  MBErrorCode result;
  MBRange proc_ents;
      // get the entities in the partition sets
  MBRange part_sets;
  result = get_partition_sets(part_sets);
  RRA(" ");
  
  for (MBRange::iterator rit = part_sets.begin(); rit != part_sets.end(); rit++) {
    MBRange tmp_ents;
    result = mbImpl->get_entities_by_handle(*rit, tmp_ents, true);
    if (MB_SUCCESS != result) return result;
    proc_ents.merge(tmp_ents);
  }

    // resolve dim is maximal dim of entities in proc_ents
  if (-1 == resolve_dim) {
    resolve_dim = mbImpl->dimension_from_handle(*proc_ents.rbegin()); 
    RRA("Couldn't get dimension.");
    
  }

    // proc_ents should all be of same dimension
  if (resolve_dim > shared_dim &&
      mbImpl->dimension_from_handle(*proc_ents.rbegin()) !=
      mbImpl->dimension_from_handle(*proc_ents.begin())) {
    MBRange::iterator lower = proc_ents.lower_bound(MBCN::TypeDimensionMap[0].first),
      upper = proc_ents.upper_bound(MBCN::TypeDimensionMap[resolve_dim-1].second);
    proc_ents.erase(lower, upper);
  }
  
    // must call even if we don't have any entities, to make sure
    // collective comm'n works
  return resolve_shared_ents(proc_ents, resolve_dim, shared_dim);
}
  
MBErrorCode MBParallelComm::resolve_shared_ents(MBRange &proc_ents,
                                                int resolve_dim,
                                                int shared_dim) 
{
  MBErrorCode result;
  if (debug) std::cerr << "Resolving shared entities." << std::endl;

  if (-1 == shared_dim) {
    if (0 == resolve_dim) {
      result = mbImpl->get_dimension(shared_dim); 
      RRA("Couldn't get dimension.");
    }
    else shared_dim = mbImpl->dimension_from_handle(*proc_ents.begin())-1;
  }
  assert(shared_dim >= 0 && resolve_dim >= 0);
  
    // get the skin entities by dimension
  MBRange skin_ents[4];
  std::vector<int> gid_data;
  std::vector<MBEntityHandle> handle_vec;
  int skin_dim;

    // get the entities to be skinned
  if (resolve_dim < shared_dim) {
      // for vertex-based partition, it's the elements adj to the vertices
    result = mbImpl->get_adjacencies(proc_ents, shared_dim,
                                     false, skin_ents[resolve_dim],
                                     MBInterface::UNION);
    RRA("Failed getting skinned entities.");
    skin_dim = shared_dim-1;
  }
  else {
      // for element-based partition, it's just the elements
    skin_ents[resolve_dim] = proc_ents;
    skin_dim = resolve_dim-1;
  }

    // find the skin
  MBSkinner skinner(mbImpl);
  result = skinner.find_skin(skin_ents[skin_dim+1], skin_ents[skin_dim],
                             skin_ents[skin_dim], true);
  RRA("Failed to find skin.");
  if (debug) std::cerr << "Found skin, now resolving." << std::endl;

    // get entities adjacent to skin ents from shared_dim down to
    // zero; don't create them if they don't exist already
  for (int this_dim = skin_dim-1; this_dim >= 0; this_dim--) {
    result = mbImpl->get_adjacencies(skin_ents[skin_dim], this_dim,
                                     false, skin_ents[this_dim],
                                     MBInterface::UNION);
    RRA("Failed getting skin adjacencies.");
  }

    // resolve shared vertices first

    // global id tag
  MBTag gid_tag; int def_val = -1;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int),
                              MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                              &def_val, true);
  if (MB_FAILURE == result) return result;

  else if (MB_ALREADY_ALLOCATED != result) {
      // just created it, so we need global ids
    result = assign_global_ids(0, skin_dim+1);
    RRA("Failed assigning global ids.");
  }

    // store index in temp tag; reuse gid_data 
  gid_data.resize(2*skin_ents[0].size());
  int idx = 0;
  for (MBRange::iterator rit = skin_ents[0].begin(); 
       rit != skin_ents[0].end(); rit++) 
    gid_data[idx] = idx, idx++;
  MBTag idx_tag;
  result = mbImpl->tag_create("__idx_tag", sizeof(int), MB_TAG_DENSE,
                              MB_TYPE_INTEGER, idx_tag, &def_val, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) return result;
  result = mbImpl->tag_set_data(idx_tag, skin_ents[0], &gid_data[0]);
  RRA("Couldn't assign index tag.");

    // get gids for skin ents in a vector, to pass to gs
  result = mbImpl->tag_get_data(gid_tag, skin_ents[0], &gid_data[0]);
  RRA("Couldn't get gid tag for skin vertices.");

    // put handles in vector for passing to gs setup
  std::copy(skin_ents[0].begin(), skin_ents[0].end(), 
            std::back_inserter(handle_vec));
  
    // get a crystal router
  crystal_data *cd = procConfig.crystal_router();
  
    // get total number of entities; will overshoot highest global id, but
    // that's ok
  int num_total, num_local;
  result = mbImpl->get_number_entities_by_dimension(0, 0, num_local);
  if (MB_SUCCESS != result) return result;
  int failure = MPI_Allreduce(&num_local, &num_total, 1,
                              MPI_INT, MPI_SUM, procConfig.proc_comm());
  if (failure) {
    result = MB_FAILURE;
    RRA("Allreduce for total number of shared ents failed.");
  }
  
    // call gather-scatter to get shared ids & procs
  gs_data *gsd;
  assert(sizeof(ulong_) == sizeof(MBEntityHandle));
  if (sizeof(int) != sizeof(ulong_)) {
    std::vector<long> lgid_data(gid_data.size());
    std::copy(gid_data.begin(), gid_data.end(), lgid_data.begin());
    gsd = gs_data_setup(skin_ents[0].size(), &lgid_data[0], 
                        (ulong_*)&handle_vec[0], 2, 1, 1, cd);
  }
  else {
    gsd = gs_data_setup(skin_ents[0].size(), (long*)&gid_data[0], 
                        (ulong_*)&handle_vec[0], 2, 1, 1, cd);
  }
  
  if (NULL == gsd) {
    result = MB_FAILURE;
    RRA("Couldn't create gs data.");
  }

    // get shared proc tags
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                sharedh_tag, sharedhs_tag, pstatus_tag);
  RRA("Couldn't get shared proc tags.");
  
    // load shared verts into a tuple, then sort by index
  tuple_list shared_verts;
  tuple_list_init_max(&shared_verts, 2, 0, 1, 0, 
                      skin_ents[0].size()*(MAX_SHARING_PROCS+1));
  unsigned int i = 0, j = 0;
  for (unsigned int p = 0; p < gsd->nlinfo->np; p++) 
    for (unsigned int np = 0; np < gsd->nlinfo->nshared[p]; np++) {
      shared_verts.vi[i++] = gsd->nlinfo->sh_ind[j];
      shared_verts.vi[i++] = gsd->nlinfo->target[p];
      shared_verts.vul[j] = gsd->nlinfo->ulabels[j];
      j++;
      shared_verts.n++;
    }
  
  int max_size = skin_ents[0].size()*(MAX_SHARING_PROCS+1);
  std::vector<int> sort_buffer(max_size);
  tuple_list_sort(&shared_verts, 0,(buffer*)&sort_buffer[0]);

    // set sharing procs and handles tags on skin ents
  int maxp = -1;
  std::vector<int> sharing_procs(MAX_SHARING_PROCS);
  std::fill(sharing_procs.begin(), sharing_procs.end(), maxp);
  j = 0; i = 0;

    // get ents shared by 1 or n procs
  std::map<std::vector<int>, MBRange> proc_nranges;
  MBRange proc_verts;
  result = mbImpl->get_adjacencies(proc_ents, 0, false, proc_verts,
                                   MBInterface::UNION);
  RRA("Couldn't get proc_verts.");
  
  result = tag_shared_verts(shared_verts, skin_ents,
                            proc_nranges, proc_verts);
  RRA("Trouble tagging shared verts.");

    // get entities shared by 1 or n procs
  result = tag_shared_ents(resolve_dim, shared_dim, shared_verts, skin_ents,
                           proc_nranges);
  RRA("Trouble tagging shared entities.");

  if (debug) {
    for (std::map<std::vector<int>, MBRange>::const_iterator mit = proc_nranges.begin();
         mit != proc_nranges.end(); mit++) {
      std::cout << "Iface: ";
      for (std::vector<int>::const_iterator vit = ((*mit).first).begin();
           vit != ((*mit).first).end(); vit++) std::cout << " " << *vit;
      std::cout << std::endl;
    }
  }
  
    // create the sets for each interface; store them as tags on
    // the interface instance
  MBRange iface_sets;
  result = create_interface_sets(proc_nranges, resolve_dim, shared_dim,
                                 &iface_sets);
  RRA("Trouble creating iface sets.");

    // resolve shared entity remote handles; implemented in ghost cell exchange
    // code because it's so similar
  result = exchange_ghost_cells(-1, -1, 0, true, true);
  RRA("Trouble resolving shared entity remote handles.");

    // now set the shared/interface tag on non-vertex entities on interface
  result = tag_iface_entities(iface_sets);
  RRA("Failed to tag iface entities.");

    // now build parent/child links for interface sets
  result = create_iface_pc_links(iface_sets);
  RRA("Trouble creating interface parent/child links.");
  
    // done
  return result;
}

MBErrorCode MBParallelComm::tag_iface_entities(MBRange &iface_ents) 
{
  MBRange all_ents, if_ents;
  MBErrorCode result;
  for (MBRange::iterator if_it = iface_ents.begin(); if_it != iface_ents.end(); if_it++) {
      // get all ents
    result = mbImpl->get_entities_by_handle(*if_it, if_ents);
    RRA("Trouble getting iface entities.");
  }
  std::vector<unsigned char> tvals(if_ents.size());
  result = mbImpl->tag_get_data(pstatus_tag(), if_ents, &tvals[0]);
  RRA("Failed to get pstatus tag values.");
  
  for (unsigned int i = 0; i < tvals.size(); i++)
    tvals[i] |= PSTATUS_INTERFACE;
  
  result = mbImpl->tag_set_data(pstatus_tag(), if_ents, &tvals[0]);
  RRA("Failed to set pstatus tag values.");
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::set_pstatus_entities(MBRange &pstatus_ents,
                                                 unsigned char pstatus_val,
                                                 bool lower_dim_ents,
                                                 bool verts_too,
                                                 int operation) 
{
  std::vector<unsigned char> pstatus_vals(pstatus_ents.size());
  MBRange all_ents, *range_ptr = &pstatus_ents;
  MBErrorCode result;
  if (lower_dim_ents || verts_too) {
    all_ents = pstatus_ents;
    range_ptr = &all_ents;
    int start_dim = (lower_dim_ents ? mbImpl->dimension_from_handle(*pstatus_ents.rbegin())-1 : 0);
    for (; start_dim >= 0; start_dim--) {
      result = mbImpl->get_adjacencies(all_ents, start_dim, true, all_ents,
                                       MBInterface::UNION);
      RRA(" ");
    }
  }
  if (MBInterface::UNION == operation) {
    result = mbImpl->tag_get_data(pstatus_tag(), *range_ptr, &pstatus_vals[0]);
    RRA("Couldn't get pstatus tag value.");
    for (unsigned int i = pstatus_ents.size()-1; i >= 0; i--)
      pstatus_vals[i] |= pstatus_val;
  }
  else {
    for (unsigned int i = pstatus_ents.size()-1; i >= 0; i--)
      pstatus_vals[i] = pstatus_val;
  }
  result = mbImpl->tag_set_data(pstatus_tag(), *range_ptr, &pstatus_vals[0]);
  RRA("Couldn't set pstatus tag value.");
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::create_interface_sets(std::map<std::vector<int>, MBRange> &proc_nranges,
                                                  int resolve_dim, int shared_dim,
                                                  MBRange *iface_sets_ptr) 
{
  if (proc_nranges.empty()) return MB_SUCCESS;
  
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag,
                                            pstatus_tag);
  RRA("Trouble getting shared proc tags in create_interface_sets.");
  MBRange::iterator rit;

  MBRange psets;
  if (!iface_sets_ptr) iface_sets_ptr = &psets;

    // get all partition sets and mark contents with iface set tag; 
    // pre-use iface_sets
  MBTag pset_tag = partition_tag();
  MBRange tmp_ents, tmp_ents2;

  result = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &pset_tag, NULL, 1,
                                                psets);
  RRA("Couldn't get PARALLEL_PARTITION sets.");

    // create interface sets, tag them, and tag their contents with iface set tag
  std::vector<MBEntityHandle> tag_vals;
  std::vector<unsigned char> pstatus;
  for (std::map<std::vector<int>,MBRange>::iterator mit = proc_nranges.begin();
       mit != proc_nranges.end(); mit++) {
      // create the set
    MBEntityHandle new_set;
    result = mbImpl->create_meshset(MESHSET_SET, new_set); 
    RRA("Failed to create interface set.");
    iface_sets_ptr->insert(new_set);

    int nump = 0;
    while (((*mit).first)[nump] != -1 && nump < MAX_SHARING_PROCS) nump++;
    
      // add entities
    result = mbImpl->add_entities(new_set, (*mit).second); 
    RRA("Failed to add entities to interface set.");
      // tag set with the proc rank(s)
    if (nump == 1)
      result = mbImpl->tag_set_data(sharedp_tag, &new_set, 1, 
                                    &((*mit).first)[0]); 
    else
      result = mbImpl->tag_set_data(sharedps_tag, &new_set, 1, 
                                    &((*mit).first)[0]); 
    RRA("Failed to tag interface set with procs.");

      // get the owning proc, then set the pstatus tag on iface set
    int min_proc = ((*mit).first)[0];
    unsigned char pval = (PSTATUS_SHARED | PSTATUS_INTERFACE);
    if (min_proc < (int) procConfig.proc_rank()) pval |= PSTATUS_NOT_OWNED;
    result = mbImpl->tag_set_data(pstatus_tag, &new_set, 1, &pval); 
    RRA("Failed to tag interface set with pstatus.");

      // tag the entities with the same thing
    pstatus.resize((*mit).second.size());
    std::fill(pstatus.begin(), pstatus.end(), pval);
    result = mbImpl->tag_set_data(pstatus_tag, (*mit).second, &pstatus[0]); 
    RRA("Failed to tag interface set entities with pstatus.");
  }

    // set tag on interface instance holding all interface sets for this instance
  MBEntityHandle tmp_iface_sets[MAX_SHARING_PROCS];
  std::fill(tmp_iface_sets, tmp_iface_sets+MAX_SHARING_PROCS, 0);
  unsigned int i;
  for (rit = iface_sets_ptr->begin(), i = 0; rit != iface_sets_ptr->end(); 
       rit++, i++)
    tmp_iface_sets[i] = *rit;
  result = mbImpl->tag_set_data(iface_sets_tag(), NULL, 0,
                                tmp_iface_sets);
  RRA("Couldn't set iface sets tag on MOAB interface.");
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::create_iface_pc_links(MBRange &iface_sets) 
{
    // now that we've resolved the entities in the iface sets, 
    // set parent/child links between the iface sets

    // first tag all entities in the iface sets
  MBTag tmp_iface_tag;
  MBEntityHandle tmp_iface_set = 0;
  MBErrorCode result = mbImpl->tag_create("__tmp_iface", sizeof(MBEntityHandle),
                                          MB_TAG_DENSE, MB_TYPE_HANDLE,
                                          tmp_iface_tag, &tmp_iface_set);
  RRA("Failed to create temporary iface set tag.");

  MBRange iface_ents;
  std::vector<MBEntityHandle> tag_vals;
  MBRange::iterator rit;
  
  for (rit = iface_sets.begin(); rit != iface_sets.end(); rit++) {
      // tag entities with interface set
    iface_ents.clear();
    result = mbImpl->get_entities_by_handle(*rit, iface_ents);
    RRA("Couldn't get entities in iface set.");
    
    if (iface_ents.empty()) continue;
    
    tag_vals.resize(iface_ents.size());
    std::fill(tag_vals.begin(), tag_vals.end(), *rit);
    result = mbImpl->tag_set_data(tmp_iface_tag, iface_ents, &tag_vals[0]); 
    RRA("Failed to tag iface entities with interface set.");
  }
  
    // now go back through interface sets and add parent/child links
  MBRange tmp_ents2;
  for (int d = 2; d >= 0; d--) {
    for (rit = iface_sets.begin(); rit != iface_sets.end(); rit++) {
        // get entities on this interface
      iface_ents.clear();
      result = mbImpl->get_entities_by_handle(*rit, iface_ents, true);
      RRA("Couldn't get entities by dimension.");
      if (iface_ents.empty() ||
          mbImpl->dimension_from_handle(*iface_ents.rbegin()) != d) continue;

        // get higher-dimensional entities and their interface sets
      result = mbImpl->get_adjacencies(&(*iface_ents.begin()), 1, d+1,
                                       false, tmp_ents2);
      RRA("Couldn't get adjacencies for interface sets.");
      tag_vals.resize(tmp_ents2.size());
      result = mbImpl->tag_get_data(tmp_iface_tag, tmp_ents2, &tag_vals[0]);
      RRA("Couldn't get iface set tag for interface sets.");
      
        // go through and for any on interface make it a parent
      MBEntityHandle last_set = 0;
      for (unsigned int i = 0; i < tag_vals.size(); i++) {
        if (tag_vals[i] && tag_vals[i] != last_set) {
          result = mbImpl->add_parent_child(tag_vals[i], *rit);
          RRA("Couldn't add parent/child link for interface set.");
          last_set = tag_vals[i];
        }
      }
    }
  }
  
    // delete the temporary tag
  result = mbImpl->tag_delete(tmp_iface_tag);
  RRA("Couldn't delete tmp iface tag.");

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::tag_shared_ents(int resolve_dim,
                                            int shared_dim,
                                            tuple_list &shared_verts,
                                            MBRange *skin_ents,
                                            std::map<std::vector<int>, MBRange> &proc_nranges) 
{
    // set sharing procs tags on other skin ents
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag, pstatus_tag);
  RRA("Trouble getting shared proc tags in tag_shared_ents.");
  const MBEntityHandle *connect; int num_connect;
  std::vector<int> sharing_procs(MAX_SHARING_PROCS);
  int sharing_procs2[MAX_SHARING_PROCS];
  std::fill(sharing_procs.begin(), sharing_procs.end(), -1);
  std::vector<unsigned char> pstatus_flags(MB_MAX_SUB_ENTITIES);

  for (int d = 3; d > 0; d--) {
    if (resolve_dim == d) continue;
    
    for (MBRange::iterator rit = skin_ents[d].begin();
         rit != skin_ents[d].end(); rit++) {
        // get connectivity
      result = mbImpl->get_connectivity(*rit, connect, num_connect);
      RRA("Failed to get connectivity on non-vertex skin entities.");
      MBRange sp_range, vp_range;

        // if any vertices not shared, this entity isn't
      assert(num_connect <= MB_MAX_SUB_ENTITIES);
      bool is_shared = true;
      result = mbImpl->tag_get_data(pstatus_tag, connect, num_connect,
                                    &pstatus_flags[0]);
      RRA("Couldn't get pstatus flag.");
      for (int nc = 0; nc < num_connect; nc++) {
        if (!(pstatus_flags[nc] & PSTATUS_SHARED)) {
          is_shared = false;
          break;
        }
      }
      if (!is_shared) continue;

      bool and_zero = false;
      for (int nc = 0; nc < num_connect; nc++) {
          // only have to initialize the 2nd one; if sharedps tag is used,
          // all get initialized
        sharing_procs2[1] = -1;
        
          // get sharing procs
        result = mbImpl->tag_get_data(sharedp_tag, connect+nc, 1, sharing_procs2);
        RRA("Couldn't get sharedp_tag on skin vertices in entity.");
        if (sharing_procs2[0] == -1) {
          result = mbImpl->tag_get_data(sharedps_tag, connect+nc, 1, sharing_procs2);
          RRA("Couldn't get sharedps_tag on skin vertices in entity.");
        }
        assert(-1 != sharing_procs2[0]);
        
          // build range of sharing procs for this vertex
        unsigned int p = 0; vp_range.clear();
        while (sharing_procs2[p] != -1 && p < MAX_SHARING_PROCS)
          vp_range.insert(sharing_procs2[p]), p++;
        assert(p < MAX_SHARING_PROCS);
          // intersect with range for this skin ent
        if (0 == nc) {
          sp_range = vp_range;
          if (sharing_procs2[0] == 0) and_zero = true;
        }
        else if (resolve_dim < shared_dim) 
          sp_range.merge(vp_range);
        else 
          sp_range = sp_range.intersect(vp_range);

          // need to also save rank zero, since ranges don't handle that
        if (and_zero && sharing_procs2[0] != 0) and_zero = false;
      }

      if (sp_range.empty() && resolve_dim < shared_dim) continue;

        // intersection is the owning proc(s) for this skin ent
      if (sp_range.empty() && !and_zero) continue;

      MBRange::iterator rit2;
        // set tag for this ent
      int j = 0;
      if (and_zero) sharing_procs[j++] = 0;
      
      for (rit2 = sp_range.begin(); rit2 != sp_range.end(); 
           rit2++, j++) sharing_procs[j] = *rit2;
      proc_nranges[sharing_procs].insert(*rit);
      if (2 > j) {
        result = mbImpl->tag_set_data(sharedp_tag, &(*rit), 1,
                                      &sharing_procs[0]);
        RRA("Failed to set sharedp_tag on non-vertex skin entity.");
      }
      else {
        result = mbImpl->tag_set_data(sharedps_tag, &(*rit), 1,
                                      &sharing_procs[0]);
        RRA("Failed to set sharedps_tag on non-vertex skin entity.");
          // set hs tag as placeholder, will get filled in later
        static std::vector<MBEntityHandle> sharing_hs(MAX_SHARING_PROCS, 0);
        result = mbImpl->tag_set_data(sharedhs_tag, &(*rit), 1,
                                      &sharing_hs[0]);
        RRA("Failed to set sharedhs_tag on non-vertex skin entity.");
      }

        // reset sharing proc(s) tags
      std::fill(sharing_procs.begin(), sharing_procs.end(), -1);
    }
  }

    // build range for each sharing proc
  std::map<int, MBRange> proc_ranges;
  for (std::map<std::vector<int>, MBRange>::iterator mit = proc_nranges.begin();
       mit != proc_nranges.end(); mit++) {
    for (unsigned int i = 0; i < (*mit).first.size(); i++) 
      proc_ranges[(*mit).first[i]].merge((*mit).second);
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::tag_shared_verts(tuple_list &shared_ents,
                                             MBRange *skin_ents,
                                             std::map<std::vector<int>, MBRange> &proc_nranges,
                                             MBRange &proc_verts) 
{
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag, pstatus_tag);
  RRA("Trouble getting shared proc tags in tag_shared_verts.");
  
  unsigned int j = 0, i = 0;
  std::vector<int> sharing_procs(MAX_SHARING_PROCS);
  MBEntityHandle sharing_handles[MAX_SHARING_PROCS];
  std::fill(sharing_procs.begin(), sharing_procs.end(), -1);
  std::fill(sharing_handles, sharing_handles+MAX_SHARING_PROCS, 0);
  
  while (j < 2*shared_ents.n) {
      // count & accumulate sharing procs
    unsigned int nump = 0;
    int this_idx = shared_ents.vi[j];
    MBEntityHandle this_ent = skin_ents[0][this_idx];
    while (shared_ents.vi[j] == this_idx) {
      j++;
      sharing_procs[nump] = shared_ents.vi[j++];
      sharing_handles[nump++] = shared_ents.vul[i++];
    }

    std::sort(&sharing_procs[0], &sharing_procs[nump]);
    proc_nranges[sharing_procs].insert(this_ent);

    if (1 == nump) {
      result = mbImpl->tag_set_data(sharedp_tag, &this_ent, 1,
                                    &sharing_procs[0]);
      result = mbImpl->tag_set_data(sharedh_tag, &this_ent, 1,
                                    sharing_handles);
    }
    else {
      result = mbImpl->tag_set_data(sharedps_tag, &this_ent, 1,
                                    &sharing_procs[0]);
      result = mbImpl->tag_set_data(sharedhs_tag, &this_ent, 1,
                                    sharing_handles);
    }
    RRA("Failed setting shared_procs tag on skin vertices.");

      // tag the entity as shared & interface too (this function only
      // ever called for interface vertices)
    unsigned char share_flag = PSTATUS_SHARED;
    result = mbImpl->tag_set_data(pstatus_tag, &this_ent, 1, &share_flag);
    RRA("Couldn't set shared tag on shared vertex.");

      // reset sharing proc(s) tags
    std::fill(sharing_procs.begin(), sharing_procs.end(), -1);
    std::fill(sharing_handles, sharing_handles+nump, 0);
  }

  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::get_iface_entities(int other_proc,
                                               int dim,
                                               MBRange &iface_ents) 
{
  MBRange iface_sets;
  std::vector<int> iface_procs;
  MBErrorCode result = get_interface_sets_procs(iface_sets, iface_procs);
  RRA("Failed to get iface sets/procs.");
  
  for (MBRange::iterator rit = iface_sets.begin(); rit != iface_sets.end(); rit++) {
    if (-1 != other_proc && !is_iface_proc(*rit, other_proc)) continue;
    
    if (-1 == dim) result = mbImpl->get_entities_by_handle(*rit, iface_ents);
    else result = mbImpl->get_entities_by_dimension(*rit, dim, iface_ents);
    RRA(" Failed to get entities in iface set.");
  }
  
  return MB_SUCCESS;
}

  //! return partition sets; if tag_name is input, gets sets with
  //! that tag name, otherwise uses PARALLEL_PARTITION tag
MBErrorCode MBParallelComm::get_partition_sets(MBRange &part_sets,
                                               const char *tag_name) 
{
  MBErrorCode result;
  
  if (NULL != tag_name) {
    result = MB_NOT_IMPLEMENTED;
    RRA("Specified tag name not yet implemented.");
  }
    
  MBTag part_tag = partition_tag();
  result = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &part_tag,
                                                NULL, 1, part_sets);
  return result;
}

  //! get communication interface sets and the processors with which
  //! this processor communicates; sets are sorted by processor
MBErrorCode MBParallelComm::get_interface_sets_procs(MBRange &iface_sets,
                                                     std::vector<int> &iface_procs)
{
  MBTag iface_tag = iface_sets_tag();
  if (0 == iface_tag) return MB_FAILURE;

    // make sure the sharing procs vector is empty
  iface_procs.clear();

    // get the iface sets, which are stored on the instance
  int tag_sz;
  MBErrorCode result = mbImpl->tag_get_size(iface_tag, tag_sz);
  RRA("Failed to get iface tag size.");
  tag_sz /= sizeof(MBEntityHandle);
  std::vector<MBEntityHandle> iface_sets_vec(tag_sz);
  result = mbImpl->tag_get_data(iface_tag, NULL, 0, &iface_sets_vec[0]);
  RRA("Failed to find iface tag.");
    // can do a straight copy because ranges don't keep 0-values
  std::copy(iface_sets_vec.begin(), iface_sets_vec.end(), 
            mb_range_inserter(iface_sets));

    // pre-load vector of single-proc tag values
  unsigned int i, j;
  std::vector<int> iface_proc(iface_sets.size());
  result = mbImpl->tag_get_data(sharedp_tag(), iface_sets, &iface_proc[0]);
  RRA("Failed to get iface_proc for iface sets.");

    // get sharing procs either from single-proc vector or by getting
    // multi-proc tag value
  int tmp_iface_procs[MAX_SHARING_PROCS];
  std::set<int> procs_set;
  std::fill(tmp_iface_procs, tmp_iface_procs+MAX_SHARING_PROCS, -1);
  MBRange::iterator rit;
  for (rit = iface_sets.begin(), i = 0; rit != iface_sets.end(); rit++, i++) {
    if (-1 != iface_proc[i]) procs_set.insert(iface_proc[i]);
    else {
        // get the sharing_procs tag
      result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1,
                                    tmp_iface_procs);
      RRA("Failed to get iface_procs for iface set.");
      for (j = 0; j < MAX_SHARING_PROCS; j++) {
        if (-1 != tmp_iface_procs[j]) procs_set.insert(tmp_iface_procs[j]);
        else {
          std::fill(tmp_iface_procs, tmp_iface_procs+j, -1);
          break;
        }
      }
    }
  }

    // now put the set contents into the vector
  std::copy(procs_set.begin(), procs_set.end(), std::back_inserter(iface_procs));
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::get_pstatus_entities(int dim,
                                                 unsigned char pstatus_val,
                                                 MBRange &pstatus_ents)
{
  MBRange ents;
  MBErrorCode result;
  
  if (-1 == dim) result = mbImpl->get_entities_by_handle(0, ents);
  else result = mbImpl->get_entities_by_dimension(0, dim, ents);
  RRA(" ");
  
  std::vector<unsigned char> pstatus(ents.size());
  result = mbImpl->tag_get_data(pstatus_tag(), ents, &pstatus[0]);
  RRA("Couldn't get pastatus tag.");
  MBRange::iterator rit = ents.begin();
  int i = 0;
  if (pstatus_val) {
    for (; rit != ents.end(); i++, rit++)
      if (pstatus[i]&pstatus_val &&
          (-1 == dim || mbImpl->dimension_from_handle(*rit) == dim)) 
        pstatus_ents.insert(*rit);
  }
  else {
    for (; rit != ents.end(); i++, rit++)
      if (!pstatus[i] &&
          (-1 == dim || mbImpl->dimension_from_handle(*rit) == dim)) 
        pstatus_ents.insert(*rit);
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::check_global_ids(MBEntityHandle this_set,
                                             const int dimension, 
                                             const int start_id,
                                             const bool largest_dim_only,
                                             const bool parallel)
{
    // global id tag
  MBTag gid_tag; int def_val = -1;
  MBErrorCode result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int),
                                          MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                                          &def_val, true);
  if (MB_ALREADY_ALLOCATED != result &&
      MB_SUCCESS != result) {
    RRA("Failed to create/get gid tag handle.");
  }

  else if (MB_ALREADY_ALLOCATED != result) {
      // just created it, so we need global ids
    result = assign_global_ids(this_set, dimension, start_id, largest_dim_only,
                               parallel);
    RRA("Failed assigning global ids.");
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_ghost_layers(MBEntityHandle iface_set,
                                             int to_dim, int bridge_dim,
                                             int num_layers,
                                             MBRange &to_ents) 
{
  MBRange bridge_ents, *bridge_ents_ptr = &bridge_ents;
  MBErrorCode result;

    // if 0 layers, put bridge_ents in to_ents, to avoid range copy later
  if (0 == num_layers) bridge_ents_ptr = &to_ents;

  if (bridge_dim == -1)
    result = mbImpl->get_entities_by_handle(iface_set, *bridge_ents_ptr);
  else
    result = mbImpl->get_entities_by_dimension(iface_set, bridge_dim, *bridge_ents_ptr);
  RRA("Couldn't get bridge ents in the set.");

  if (bridge_ents_ptr->empty()) return MB_SUCCESS;


    // for each layer, get bridge-adj entities and accumulate
  for (int nl = 0; nl < num_layers; nl++) {
    MBRange new_ghosts, new_bridges;
    result = mbImpl->get_adjacencies(*bridge_ents_ptr, to_dim, false, new_ghosts,
                                     MBInterface::UNION);
    RRA("Trouble getting ghost adjacencies in iteration.");
    to_ents.merge(new_ghosts);
    result = mbImpl->get_adjacencies(new_ghosts, bridge_dim, false, new_bridges,
                                     MBInterface::UNION);
    if (nl < num_layers-1)
      *bridge_ents_ptr = new_bridges.subtract(*bridge_ents_ptr);
  }

  return MB_SUCCESS;
}

bool MBParallelComm::is_iface_proc(MBEntityHandle this_set,
                                   int to_proc) 
{
  int sharing_procs[MAX_SHARING_PROCS];
  std::fill(sharing_procs, sharing_procs+MAX_SHARING_PROCS, -1);
  MBErrorCode result = mbImpl->tag_get_data(sharedp_tag(), &this_set, 1,
                                            sharing_procs);
  if (to_proc == sharing_procs[0]) return true;
  
  result = mbImpl->tag_get_data(sharedps_tag(), &this_set, 1,
                                sharing_procs);
  for (int i = 0; i < MAX_SHARING_PROCS; i++) {
    if (to_proc == sharing_procs[i]) return true;
    else if (-1 == sharing_procs[i]) return false;
  }
  
  return false;
}

MBErrorCode MBParallelComm::remove_nonowned_shared(MBRange &ents,
                                                   int to_proc,
                                                   bool owned_test) 
{
    // remove from ents any entities which are not owned locally or
    // who are already shared with to_proc
  std::vector<unsigned char> shared_flags(ents.size());
  MBErrorCode result = mbImpl->tag_get_data(pstatus_tag(), ents,
                                            &shared_flags[0]);
  RRA("Failed to get pstatus flag.");
  MBRange tmp_ents;
  MBRange::iterator rit;
  int sharing_procs[MAX_SHARING_PROCS];
  std::fill(sharing_procs, sharing_procs+MAX_SHARING_PROCS, -1);
  int i;
  for (rit = ents.begin(), i = 0; rit != ents.end(); rit++, i++) {
      // don't save if I don't own this entity
    if (owned_test && (PSTATUS_NOT_OWNED & shared_flags[i])) continue;
  
      // if entity isn't shared yet, add to list
    else if (!(PSTATUS_SHARED & shared_flags[i])) tmp_ents.insert(*rit);
    
      // else we need to check sharing procs
    else {
      result = mbImpl->tag_get_data(sharedp_tag(), &(*rit), 1,
                                    sharing_procs);
      RRA(" ");
      if (-1 == sharing_procs[0]) {
        result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1,
                                      sharing_procs);
        assert(-1 != sharing_procs[0]);
        RRA(" ");
      }
      for (unsigned int j = 0; j < MAX_SHARING_PROCS; j++) {
          // if to_proc shares this entity, skip it
        if (sharing_procs[j] == to_proc) break;
          // if we get here, no more sharing procs, and it's not shared
          // with to_proc, so add to list
        else if (sharing_procs[j] == -1) 
          tmp_ents.insert(*rit);
      }
      std::fill(sharing_procs, sharing_procs+MAX_SHARING_PROCS, -1);
    }
  }
    
  ents.swap(tmp_ents);
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::exchange_ghost_cells(int ghost_dim, int bridge_dim,
                                                 int num_layers,
                                                 bool store_remote_handles,
                                                 bool wait_all)
{
    // if we're only finding out about existing ents, we have to be storing
    // remote handles too
  assert(num_layers > 0 || store_remote_handles);
  
    // get the b-dimensional interface(s) with with_proc, where b = bridge_dim
  MBErrorCode result;
  int success;

    // get all procs interfacing to this proc
  MBRange iface_sets;
  std::vector<int> iface_procs;
  result = get_interface_sets_procs(iface_sets, iface_procs);
  RRA("Failed to get iface sets, procs");

    // post ghost irecv's for all interface procs
    // index greqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_reqs(2*MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  std::vector<MPI_Status> gstatus(MAX_SHARING_PROCS);
  std::vector<int>::iterator vit;
  for (vit = iface_procs.begin(); vit != iface_procs.end(); vit++) {
    int ind = get_buffers(*vit);
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, *vit,
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    // pack and send ghosts from this proc to others
    // make sendReqs vector to simplify initialization
  std::fill(sendReqs, sendReqs+2*MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  MBRange recd_ents[MAX_SHARING_PROCS], sent_ents[MAX_SHARING_PROCS];
  
    // keep track of new ghosted and ghost ents so we can tag them later
  MBRange new_ghosted, new_ghosts;
  
  for (vit = iface_procs.begin(); vit != iface_procs.end(); vit++) {
    int ind = get_buffers(*vit);
    
    MBRange bridge_ents;

      // get bridge ents on interface(s)
    for (MBRange::iterator rit = iface_sets.begin(); rit != iface_sets.end();
         rit++) {
      if (!is_iface_proc(*rit, *vit)) continue;
      
      result = get_ghost_layers(*rit, ghost_dim, bridge_dim, num_layers, bridge_ents);
      RRA("Failed to get ghost layers.");
    }

    if (0 == num_layers) {
        // just sending iface ents - take out vertices, since these are always known
      std::pair<MBRange::const_iterator, MBRange::const_iterator> its =
          bridge_ents.equal_range(MBVERTEX);
      bridge_ents.erase((MBRange::iterator)its.first,
                        (MBRange::iterator)its.second);
    }

      // pack-send; this also posts receives if store_remote_handles is true
    result = pack_send_entities(*vit, bridge_ents, false, true, 
                                store_remote_handles, (0 == num_layers),
                                ownerSBuffs[ind], ownerRBuffs[MAX_SHARING_PROCS+ind], 
                                sendReqs[ind], recv_reqs[MAX_SHARING_PROCS+ind], 
                                sent_ents[ind]);
    RRA("Failed to pack-send in ghost exchange.");

    if (0 != num_layers) new_ghosted.merge(sent_ents[ind]);
  }
  
    // receive/unpack entities
    // number of incoming messages depends on whether we're getting back
    // remote handles
  int num_incoming = iface_sets.size() * (store_remote_handles ? 2 : 1);
  
  while (num_incoming) {
    int ind;
    MPI_Status status;
    success = MPI_Waitany(2*MAX_SHARING_PROCS, &recv_reqs[0], &ind, &status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
    std::vector<MBEntityHandle> remote_handles_v, sent_ents_tmp;
    MBRange remote_handles_r;
    int new_size, num_eh;
    unsigned char *buff_ptr;
    
      // branch on message type
    switch (status.MPI_TAG) {
      case MB_MESG_SIZE:
          // incoming message just has size; resize buffer and re-call recv,
          // then re-increment incoming count
        assert(ind < MAX_SHARING_PROCS);
        new_size = *((int*)&ghostRBuffs[ind][0]);
        assert(0 > new_size);
        result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind]);
        RRA("Failed to resize recv buffer.");
        num_incoming++;
        break;
      case MB_MESG_ENTS:
          // incoming ghost entities; process
        result = recv_unpack_entities(buffProcs[ind], store_remote_handles,
                                      (0 == num_layers), 
                                      ghostRBuffs[ind], ghostSBuffs[ind], 
                                      sendReqs[ind], recd_ents[ind]);
        RRA("Failed to recv-unpack message.");
        if (0 != num_layers) new_ghosts.merge(recd_ents[ind]);
        break;
      case MB_MESG_REMOTE_HANDLES_VECTOR:
          // incoming remote handles; use to set remote handles
          // this version only called if we're exchanging iface layer...
        assert(0 == num_layers);
        buff_ptr = &ownerRBuffs[ind][0];
        assert(ind >= MAX_SHARING_PROCS);
        UNPACK_INT(buff_ptr, num_eh);
        remote_handles_v.resize(num_eh);
        UNPACK_EH(buff_ptr, &remote_handles_v[0], num_eh);
        std::copy(sent_ents[ind-MAX_SHARING_PROCS].begin(),
                  sent_ents[ind-MAX_SHARING_PROCS].end(),
                  std::back_inserter(sent_ents_tmp));
        assert(!sent_ents_tmp.empty());
        result = set_remote_data(&sent_ents_tmp[0], &remote_handles_v[0], sent_ents_tmp.size(),
                                 buffProcs[ind-MAX_SHARING_PROCS]);
        RRA("Trouble setting remote data range on sent entities in ghost exchange.");
        result = update_iface_sets(iface_sets, sent_ents[ind-MAX_SHARING_PROCS],
                                   remote_handles_v, buffProcs[ind-MAX_SHARING_PROCS]);
        RRA("Trouble updating iface sets.");
        break;
      case MB_MESG_REMOTE_HANDLES_RANGE:
          // incoming remote handles; use to set remote handles
          // this version only called if we're NOT exchanging iface layer...
        assert(0 != num_layers);
        buff_ptr = &ownerRBuffs[ind][0];
        assert(ind >= MAX_SHARING_PROCS);
        UNPACK_RANGE(buff_ptr, remote_handles_r);
        assert(!sent_ents[ind-MAX_SHARING_PROCS].empty());
        result = set_remote_data(sent_ents[ind-MAX_SHARING_PROCS], remote_handles_r, 
                                 buffProcs[ind-MAX_SHARING_PROCS]);
        RRA("Trouble setting remote data range on sent entities in ghost exchange.");
        break;
    }
  }
  
  if (0 != num_layers) {
    unsigned char pval = (PSTATUS_SHARED | PSTATUS_GHOST | PSTATUS_NOT_OWNED);
    std::vector<unsigned char> pvals(MAX(new_ghosts.size(), new_ghosted.size()), pval);
      // new ghost entities are always new
    result = mbImpl->tag_set_data(pstatus_tag(), new_ghosts, &pvals[0]);
    RRA("Trouble setting tags for new ghost entities");
    
    pval = PSTATUS_SHARED;
    std::fill(pvals.begin(), pvals.end(), pval);
    result = mbImpl->tag_set_data(pstatus_tag(), new_ghosted, &pvals[0]);
    RRA("Trouble setting tags for new ghosted entities");
  }
  
    // ok, now wait if requested
  if (wait_all) {
    MPI_Status status[2*MAX_SHARING_PROCS];
    success = MPI_Waitall(2*MAX_SHARING_PROCS, &sendReqs[0], status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failure in waitall in ghost exchange.");
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::update_iface_sets(MBRange &iface_sets, 
                                              MBRange &sent_ents,
                                              std::vector<MBEntityHandle> &remote_handles, 
                                              int from_proc) 
{
  std::vector<MBEntityHandle>::iterator remote_it = remote_handles.begin();
  MBRange::iterator sent_it = sent_ents.begin();
  MBRange ents_to_remove;
  for (; sent_it != sent_ents.end(); sent_it++, remote_it++) {
    if (!*remote_it) ents_to_remove.insert(*sent_it);
  }
  
  for (MBRange::iterator set_it = iface_sets.begin(); set_it != iface_sets.end(); set_it++) {
    if (!is_iface_proc(*set_it, from_proc)) continue;
    MBErrorCode result = mbImpl->remove_entities(*set_it, ents_to_remove);
    RRA("Couldn't remove entities from iface set in update_iface_sets.");
  }
  
  return MB_SUCCESS;
}

  //! return sharedp tag
MBTag MBParallelComm::sharedp_tag()
{
  if (!sharedpTag) {
    int def_val = -1;
    MBErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_PROC_TAG_NAME, 
                                            sizeof(int), 
                                            MB_TAG_DENSE,
                                            MB_TYPE_INTEGER, sharedpTag, 
                                            &def_val, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
      return 0;
  }
  
  return sharedpTag;
}

  //! return sharedps tag
MBTag MBParallelComm::sharedps_tag()
{
  if (!sharedpsTag) {
    MBErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_PROCS_TAG_NAME, 
                                            MAX_SHARING_PROCS*sizeof(int), 
                                            MB_TAG_SPARSE,
                                            MB_TYPE_INTEGER, sharedpsTag, NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
      return 0;
  }
  
  return sharedpsTag;
}
  
  //! return sharedh tag
MBTag MBParallelComm::sharedh_tag()
{
  if (!sharedhTag) {
    int def_val = 0;
    MBErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_HANDLE_TAG_NAME, 
                                            sizeof(MBEntityHandle), 
                                            MB_TAG_DENSE,
                                            MB_TYPE_INTEGER, sharedhTag, 
                                            &def_val, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return sharedhTag;
}
  
  //! return sharedhs tag
MBTag MBParallelComm::sharedhs_tag()
{  
  if (!sharedhsTag) {
    MBErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_HANDLES_TAG_NAME, 
                                            MAX_SHARING_PROCS*sizeof(MBEntityHandle), 
                                            MB_TAG_SPARSE,
                                            MB_TYPE_INTEGER, sharedhsTag, NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
      return 0;
  }

  return sharedhsTag;
}
  
  //! return pstatus tag
MBTag MBParallelComm::pstatus_tag()
{  
  if (!pstatusTag) {
    unsigned char tmp_pstatus = 0;
    MBErrorCode result = mbImpl->tag_create(PARALLEL_STATUS_TAG_NAME, 
                                            sizeof(unsigned char),
                                            MB_TAG_DENSE,
                                            MB_TYPE_OPAQUE, pstatusTag, 
                                            &tmp_pstatus, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return pstatusTag;
}
  
  //! return iface set tag
MBTag MBParallelComm::iface_sets_tag()
{  
  if (!ifaceSetsTag) {
    MBErrorCode result = mbImpl->tag_create(PARALLEL_IFACE_SETS_TAG_NAME, 
                                            MAX_SHARING_PROCS*sizeof(MBEntityHandle),
                                            MB_TAG_SPARSE,
                                            MB_TYPE_HANDLE, ifaceSetsTag, 
                                            NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return ifaceSetsTag;
}
  
  //! return partition set tag
MBTag MBParallelComm::partition_tag()
{  
  if (!partitionTag) {
    MBErrorCode result = mbImpl->tag_create(PARALLEL_PARTITION_TAG_NAME, 
                                            MAX_SHARING_PROCS*sizeof(MBEntityHandle),
                                            MB_TAG_SPARSE,
                                            MB_TYPE_HANDLE, partitionTag, 
                                            NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return partitionTag;
}
  
#ifdef TEST_PARALLELCOMM

#include <iostream>

#include "MBCore.hpp"
#include "MBParallelComm.hpp"
#include "MBRange.hpp"

#define PM {std::cerr << "Test failed; error message:" << std::endl;\
          std::string errmsg; \
          dynamic_cast<MBCore*>(my_impl)->get_last_error(errmsg); \
          std::cerr << errmsg << std::endl;\
          return 1;}

int main(int argc, char* argv[])
{

    // Check command line arg
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " <mesh_file_name>" << std::endl;
    exit(1);
  }

  const char* file = argv[1];
  MBCore *my_impl = new MBCore(0, 2);
  MBInterface* mbImpl = my_impl;

    // create a communicator class, which will start mpi too
  MBParallelComm pcomm(mbImpl, my_impl->tag_server(), my_impl->sequence_manager());
  MBErrorCode result;

    // load the mesh
  result = mbImpl->load_mesh(file, 0, 0);
  if (MB_SUCCESS != result) return result;

    // get the mesh
  MBRange all_mesh, whole_range;
  result = mbImpl->get_entities_by_dimension(0, 3, all_mesh);
  if (MB_SUCCESS != result) return result;
    
  int buff_size;
  result = pcomm.pack_buffer(all_mesh, false, true, true, false, whole_range, buff_size);
  PM;


    // allocate space in the buffer
  pcomm.buffer_size(buff_size);

    // pack the actual buffer
  int actual_buff_size;
  result = pcomm.pack_buffer(whole_range, false, true, false, false, all_mesh, 
                             actual_buff_size);
  PM;

    // list the entities that got packed
  std::cout << "ENTITIES PACKED:" << std::endl;
  mbImpl->list_entities(all_mesh);

    // get the buffer
  std::vector<unsigned char> tmp_buffer;
  pcomm.take_buffer(tmp_buffer);
    
    // stop and restart MOAB
  delete mbImpl;
  my_impl = new MBCore(1, 2);
  mbImpl = my_impl;
    
    // create a new communicator class, using our old buffer
  MBParallelComm pcomm2(mbImpl, my_impl->tag_server(), my_impl->sequence_manager(),
                        tmp_buffer);

    // unpack the results
  all_mesh.clear();
  result = pcomm2.unpack_buffer(all_mesh, store_remote_handles, from_proc);
  PM;
  
  std::cout << "ENTITIES UNPACKED:" << std::endl;
  mbImpl->list_entities(all_mesh);
  
  std::cout << "Success, processor " << mbImpl->proc_rank() << "." << std::endl;
  
  return 1;
}
#endif
