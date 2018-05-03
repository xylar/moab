/*
 * ParCommGraph.cpp
 *
 */

#include "moab/ParCommGraph.hpp"
// we need to recompute adjacencies for merging to work
#include "moab/Core.hpp"
#include "AEntityFactory.hpp"
namespace moab {

ParCommGraph::ParCommGraph(MPI_Comm joincomm, MPI_Group group1, MPI_Group group2, int coid1, int coid2):
comm(joincomm), compid1(coid1), compid2(coid2) {
  // find out the tasks from each group, in the joint communicator
  find_group_ranks(group1, comm, senderTasks);
  find_group_ranks(group2, comm, receiverTasks);

  rootSender = rootReceiver = false;
  rankInGroup1 =  rankInGroup2 = rankInJoin = -1; // not initialized, or not part of the group
  int mpierr = MPI_Group_rank(group1, &rankInGroup1);
  if (MPI_SUCCESS != mpierr || rankInGroup1 == MPI_UNDEFINED)
    rankInGroup1 = -1;
  mpierr = MPI_Group_rank(group2, &rankInGroup2);
  if (MPI_SUCCESS != mpierr || rankInGroup2 == MPI_UNDEFINED)
    rankInGroup2 = -1;
  mpierr = MPI_Comm_rank(comm, &rankInJoin);
  if (MPI_SUCCESS != mpierr) // it should be a fatal error
    rankInJoin = -1;
  mpierr = MPI_Comm_size(comm, &joinSize);
  if (MPI_SUCCESS != mpierr) // it should be a fatal error
    joinSize = -1;

  if (0==rankInGroup1)rootSender=true;
  if (0==rankInGroup2)rootReceiver=true;

  comm_graph = NULL;
}

ParCommGraph::~ParCommGraph() {
  // TODO Auto-generated destructor stub
}

// utility to find out the ranks of the processes of a group, with respect to a joint comm,
// which spans for sure the group
// it is used locally (in the constructor), but it can be used as a utility
void ParCommGraph::find_group_ranks(MPI_Group group, MPI_Comm joincomm, std::vector<int> & ranks)
{
   MPI_Group global_grp;
   MPI_Comm_group(joincomm, &global_grp);

   int grp_size;

   MPI_Group_size(group, &grp_size);
   std::vector<int> rks (grp_size);
   ranks.resize(grp_size);

   for (int i = 0; i < grp_size; i++)
     rks[i] = i;

   MPI_Group_translate_ranks(group, grp_size, &rks[0], global_grp, &ranks[0]);
   MPI_Group_free(&global_grp);
   return;
}

ErrorCode ParCommGraph::compute_trivial_partition (std::vector<int> & numElemsPerTaskInGroup1)
{

  recv_graph.clear(); recv_sizes.clear();
  sender_graph.clear(); sender_sizes.clear();

  if (numElemsPerTaskInGroup1.size() != senderTasks.size())
    return MB_FAILURE; // each sender has a number of elements that it owns

  // first find out total number of elements to be sent from all senders
  int total_elems=0;
  std::vector<int> accum;
  accum.push_back(0);

  int num_senders = (int) senderTasks.size();

  for (size_t k=0; k<numElemsPerTaskInGroup1.size(); k++)
  {
    total_elems+=numElemsPerTaskInGroup1[k];
    accum.push_back(total_elems);
  }

  int num_recv  =  ((int)receiverTasks.size());
  // in trivial partition, every receiver should get about total_elems/num_receivers elements
  int num_per_receiver = (int)(total_elems/num_recv);
  int leftover = total_elems - num_per_receiver*num_recv;

  // so receiver k will receive  [starts[k], starts[k+1] ) interval
  std::vector<int> starts;
  starts.resize(num_recv+1);
  starts[0]=0;
  for (int k=0; k<num_recv; k++)
  {
    starts[k+1] = starts[k]+  num_per_receiver;
    if (k<leftover) starts[k+1]++;
  }

  // each sender will send to a number of receivers, based on how the
  // arrays starts[0:num_recv] and accum[0:sendr] overlap
  int lastUsedReceiverRank = 0; // first receiver was not treated yet
  for (int j = 0; j < num_senders; j++ )
  {
    // we could start the receiver loop with the latest receiver that received from previous sender
    for (int k = lastUsedReceiverRank; k<num_recv; k++ )
    {
      // if overlap:
      if (starts[k]<accum[j+1] && starts[k+1]> accum[j] )
      {
        recv_graph[receiverTasks[k]].push_back(senderTasks[j]);
        sender_graph[senderTasks[j]].push_back(receiverTasks[k]);

        // we still need to decide what is the overlap
        int sizeOverlap = 1; // at least 1, for sure
        //1
        if ( starts[k] >= accum[j]) // one end is starts[k]
        {
          if (starts[k+1] >= accum[j+1]) // the other end is accum[j+1]
            sizeOverlap = accum[j+1]-starts[k];
          else //
            sizeOverlap = starts[k+1]-starts[k];
        }
        else // one end is accum[j]
        {
          if (starts[k+1] >= accum[j+1]) // the other end is accum[j+1]
            sizeOverlap = accum[j+1]-accum[j];
          else
            sizeOverlap = starts[k+1]-accum[j];
        }
        recv_sizes[receiverTasks[k]].push_back(sizeOverlap); // basically, task k will receive from
                                                        //   sender j, sizeOverlap elems
        sender_sizes[senderTasks[j]].push_back(sizeOverlap);
        if (starts[k] > accum[j+1])
        {
          lastUsedReceiverRank = k-1; // so next k loop will start a little higher, we probably
                                            // finished with first few receivers (up to receiver lastUsedReceiverRank)
          break ; // break the k loop, we distributed all elements from sender j to some receivers
        }
      }
    }
  }

  return MB_SUCCESS;
}
ErrorCode ParCommGraph::pack_receivers_graph(std::vector<int> & packed_recv_array )
{
  // it will basically look at local data, to pack communication graph, each receiver task will have to post receives
  // for each sender task that will send data to it;
  // the array will be communicated to root receiver, and eventually distributed to receiver tasks

  /*
   * packed_array will have receiver, number of senders, then senders, etc
   */
  for (std::map<int, std::vector<int> >::iterator it=recv_graph.begin(); it!=recv_graph.end(); it++ )
  {
    int recv = it->first;
    std::vector<int> & senders = it->second;
    packed_recv_array.push_back(recv);
    packed_recv_array.push_back( (int) senders.size() );

    for (int k = 0; k<(int)senders.size(); k++ )
      packed_recv_array.push_back( senders[k] );
  }

  return MB_SUCCESS;
}

ErrorCode ParCommGraph::split_owned_range (int sender_rank, Range & owned, std::map<int, Range> & split_ranges)
{
  int senderTask = senderTasks[sender_rank];
  std::vector<int> & distribution = sender_sizes[senderTask];
  std::vector<int> & receivers = sender_graph[senderTask];
  if (distribution.size() != receivers.size()) //
    return MB_FAILURE;

  Range current = owned; // get the full range first, then we will subtract stuff, for
  // the following ranges

  Range rleftover=current;
  for (size_t k=0; k<receivers.size(); k++ )
  {
    Range newr;
    newr.insert(current.begin(), current.begin() +distribution[k]);
    split_ranges[ receivers[k] ] = newr;

    rleftover = subtract(current, newr );
    current = rleftover;
  }

  return MB_SUCCESS;
}

// use for this the corresponding tasks and sizes
ErrorCode ParCommGraph::split_owned_range_general (Range & owned, std::map<int, Range> & split_ranges)
{
  if (corr_tasks.size() != corr_sizes.size()) //
    return MB_FAILURE;

  Range current = owned; // get the full range first, then we will subtract stuff, for
  // the following ranges

  Range rleftover=current;
  for (size_t k=0; k<corr_tasks.size(); k++ )
  {
    Range newr;
    newr.insert(current.begin(), current.begin() +corr_sizes[k]);
    split_ranges[ corr_tasks[k] ] = newr;

    rleftover = subtract(current, newr );
    current = rleftover;
  }

  return MB_SUCCESS;
}

ErrorCode ParCommGraph::send_graph(MPI_Comm jcomm)
{
  if ( is_root_sender())
  {
    int ierr;
    // will need to build a communication graph, because each sender knows now to which receiver to send data
    // the receivers need to post receives for each sender that will send data to them
    // will need to gather on rank 0 on the sender comm, global ranks of sender with receivers to send
    // build communication matrix, each receiver will receive from what sender

    std::vector<int> packed_recv_array;
    ErrorCode rval = pack_receivers_graph(packed_recv_array );
    if (MB_SUCCESS!=rval) return rval;

    int size_pack_array = (int) packed_recv_array.size();
    comm_graph = new int[size_pack_array+1];
    comm_graph[0]=size_pack_array;
    for (int k=0; k<size_pack_array; k++)
      comm_graph[k+1]=packed_recv_array[k];
    // will add 2 requests
    /// use tag 10 to send size and tag 20 to send the packed array
    // will allocate some send requests
    int nsendsroot = 2*(int) sender_graph[senderTasks[0]].size() + 2; // also need to send the graph
    sendReqs.resize(nsendsroot);
    ierr = MPI_Isend(&comm_graph[0], 1, MPI_INT, receiver(0), 10, jcomm, &sendReqs[0]); // we have to use global communicator
    if (ierr!=0) return MB_FAILURE;
    ierr = MPI_Isend(&comm_graph[1], size_pack_array, MPI_INT, receiver(0), 20, jcomm, &sendReqs[1]); // we have to use global communicator
    if (ierr!=0) return MB_FAILURE;
  }
  return MB_SUCCESS;
}

// pco has MOAB too get_moab()
ErrorCode ParCommGraph::send_mesh_parts(MPI_Comm jcomm, ParallelComm * pco, Range & owned )
{
  std::map<int, Range> split_ranges;
  ErrorCode rval = split_owned_range (rankInGroup1, owned, split_ranges);
  if (rval!=MB_SUCCESS) return rval;

  // we know this on the sender side:
  corr_tasks = sender_graph[ senderTasks[rankInGroup1]]; // copy
  corr_sizes = sender_sizes[ senderTasks[rankInGroup1]]; // another copy
  int indexReq=0;
  int ierr; // MPI error
  if (is_root_sender()) indexReq = 2; // for sendReqs
  sendReqs.resize(indexReq+2*split_ranges.size());
  for (std::map<int, Range>::iterator it=split_ranges.begin(); it!=split_ranges.end(); it++)
  {
    int receiver_proc = it->first;
    Range ents = it->second;

    // add necessary vertices too
    Range verts;
    rval = pco->get_moab()->get_adjacencies(ents, 0, false, verts, Interface::UNION);
    if (rval!=MB_SUCCESS) return rval;
    ents.merge(verts);
    ParallelComm::Buffer * buffer = new ParallelComm::Buffer(ParallelComm::INITIAL_BUFF_SIZE);
    buffer->reset_ptr(sizeof(int));
    rval = pco->pack_buffer(ents, false, true, false, -1, buffer);
    if (rval!=MB_SUCCESS) return rval;
    int size_pack = buffer->get_current_size();

    // TODO there could be an issue with endian things; check !!!!!
    // we are sending the size of the buffer first as an int!!!
    ierr = MPI_Isend(buffer->mem_ptr, 1, MPI_INT, receiver_proc, 1, jcomm, &sendReqs[indexReq]); // we have to use global communicator
    if (ierr!=0) return MB_FAILURE;
    indexReq++;

    ierr = MPI_Isend(buffer->mem_ptr, size_pack, MPI_CHAR, receiver_proc, 2, jcomm, &sendReqs[indexReq]); // we have to use global communicator
    if (ierr!=0) return MB_FAILURE;
    indexReq++;
    localSendBuffs.push_back(buffer);

  }
  return MB_SUCCESS;
}
// this is called on receiver side
ErrorCode ParCommGraph::receive_comm_graph(MPI_Comm jcomm, ParallelComm *pco, std::vector<int> & pack_array)
{
  // first, receive from sender_rank 0, the communication graph (matrix), so each receiver
  // knows what data to expect
  MPI_Comm receive = pco->comm();
  int size_pack_array, ierr;
  MPI_Status status;
  if (rootReceiver)
  {
    ierr = MPI_Recv (&size_pack_array, 1, MPI_INT, sender(0), 10, jcomm, &status);
    if (0!=ierr) return MB_FAILURE;
#ifdef VERBOSE
    std::cout <<" receive comm graph size: " << size_pack_array << "\n";
#endif
    pack_array.resize (size_pack_array);
    ierr = MPI_Recv (&pack_array[0], size_pack_array, MPI_INT, sender(0), 20, jcomm, &status);
    if (0!=ierr) return MB_FAILURE;
#ifdef VERBOSE
    std::cout <<" receive comm graph " ;
    for (int k=0; k<(int)pack_array.size(); k++)
      std::cout << " " << pack_array[k];
    std::cout <<"\n";
#endif
  }

  // now broadcast this whole array to all receivers, so they know what to expect
  ierr = MPI_Bcast(&size_pack_array, 1, MPI_INT, 0, receive);
  if (0!=ierr) return MB_FAILURE;
  pack_array.resize(size_pack_array);
  ierr = MPI_Bcast(&pack_array[0], size_pack_array, MPI_INT, 0, receive);
  if (0!=ierr) return MB_FAILURE;
  return MB_SUCCESS;
}
ErrorCode ParCommGraph::receive_mesh(MPI_Comm jcomm, ParallelComm *pco, EntityHandle local_set, std::vector<int> &senders_local)
{
  ErrorCode rval;
  int ierr;
  MPI_Status status;
  // we also need to fill corresponding mesh info on the other side
  corr_tasks = senders_local;
  Range newEnts;

  Tag orgSendProcTag; // this will be a tag set on the received mesh, with info about from what task / PE the
  // primary element came from, in the joint communicator ; this will be forwarded by coverage mesh
  int defaultInt=-1; // no processor, so it was not migrated from somewhere else
  rval = pco->get_moab()->tag_get_handle("orig_sending_processor", 1, MB_TYPE_INTEGER, orgSendProcTag,
      MB_TAG_DENSE | MB_TAG_CREAT, &defaultInt);MB_CHK_SET_ERR(rval, "can't create original sending processor tag");
  if (!senders_local.empty())
  {
    for (size_t k=0; k< senders_local.size(); k++)
    {
      int sender1 = senders_local[k]; // first receive the size of the buffer

      int size_pack;
      ierr = MPI_Recv (&size_pack, 1, MPI_INT, sender1, 1, jcomm, &status);
      if (0!=ierr) return MB_FAILURE;
      // now resize the buffer, then receive it
      ParallelComm::Buffer * buffer = new ParallelComm::Buffer(size_pack);
      //buffer->reserve(size_pack);

      ierr = MPI_Recv (buffer->mem_ptr, size_pack, MPI_CHAR, sender1, 2, jcomm, &status);
      if (0!=ierr) return MB_FAILURE;
      // now unpack the buffer we just received
      Range entities;
      std::vector<std::vector<EntityHandle> > L1hloc, L1hrem;
      std::vector<std::vector<int> > L1p;
      std::vector<EntityHandle> L2hloc, L2hrem;
      std::vector<unsigned int> L2p;

      buffer->reset_ptr(sizeof(int));
      std::vector<EntityHandle> entities_vec(entities.size());
      std::copy(entities.begin(), entities.end(), entities_vec.begin());
      rval = pco->unpack_buffer(buffer->buff_ptr, false, -1, -1, L1hloc, L1hrem, L1p, L2hloc,
                                  L2hrem, L2p, entities_vec);
      delete buffer;
      if (MB_SUCCESS!= rval) return rval;

      std::copy(entities_vec.begin(), entities_vec.end(), range_inserter(entities));
      // we have to add them to the local set
      rval = pco->get_moab()->add_entities(local_set, entities);
      if (MB_SUCCESS!= rval) return rval;
      // corr_sizes is the size of primary entities received
      Range verts = entities.subset_by_dimension(0);
      Range local_primary_ents=subtract(entities, verts);
      corr_sizes.push_back( (int) local_primary_ents.size());
      // set a tag with the original sender for the primary entity
      // will be used later for coverage mesh
      std::vector<int> orig_senders(local_primary_ents.size(), sender1);
      rval = pco->get_moab()->tag_set_data(orgSendProcTag, local_primary_ents, &orig_senders[0]);
      newEnts.merge(entities);

#ifdef VERBOSE
      std::ostringstream partial_outFile;

      partial_outFile <<"part_send_" <<sender1<<"."<< "recv"<< rankInJoin <<".vtk";

        // the mesh contains ghosts too, but they are not part of mat/neumann set
        // write in serial the file, to see what tags are missing
      std::cout<< " writing from receiver " << rankInJoin << " from sender " <<  sender1 << " entities: " << entities.size() << std::endl;
      rval = pco->get_moab()->write_file(partial_outFile.str().c_str(), 0, 0, &local_set, 1); // everything on local set received
      if (MB_SUCCESS!= rval) return rval;
#endif
    }

  }
  // make sure adjacencies are updated on the new elements

  // in order for the merging to work, we need to be sure that the adjacencies are updated (created)
  Range local_verts = newEnts.subset_by_type(MBVERTEX);
  newEnts = subtract(newEnts, local_verts);
  Core * mb = (Core*)pco->get_moab();
  AEntityFactory* adj_fact = mb->a_entity_factory();
  if (!adj_fact->vert_elem_adjacencies())
    adj_fact->create_vert_elem_adjacencies();
  else
  {
    for (Range::iterator it=newEnts.begin(); it!=newEnts.end(); it++)
    {
      EntityHandle eh = *it;
      const EntityHandle * conn = NULL;
      int num_nodes=0;
      rval = mb->get_connectivity(eh, conn, num_nodes);
      if (MB_SUCCESS!= rval) return rval;
      adj_fact->notify_create_entity(eh, conn, num_nodes);
    }
  }

  return MB_SUCCESS;
}
ErrorCode ParCommGraph::release_send_buffers(MPI_Comm jcomm)
{

  int ierr, nsize = (int)sendReqs.size();
  std::vector<MPI_Status> mult_status;
  mult_status.resize(sendReqs.size());
  ierr = MPI_Waitall(nsize, &sendReqs[0], &mult_status[0]);

  if (ierr!=0)
    return MB_FAILURE;
  // now we can free all buffers
  delete [] comm_graph;
  comm_graph = NULL;
  std::vector<ParallelComm::Buffer*>::iterator vit;
  for (vit = localSendBuffs.begin(); vit != localSendBuffs.end(); ++vit)
    delete (*vit);
  localSendBuffs.clear();
  return MB_SUCCESS;
}
// again, will use the send buffers, for nonblocking sends;
// should be the receives non-blocking too?
ErrorCode ParCommGraph::send_tag_values (MPI_Comm jcomm, ParallelComm *pco, Range & owned, Tag & tag_handle )
{
  // basically, owned.size() needs to be equal to sum(corr_sizes)
  // get info about the tag size, type, etc
  int ierr;
  Core * mb = (Core*)pco->get_moab();
  // get info about the tag
  //! Get the size of the specified tag in bytes
  int bytes_per_tag; // we need to know, to allocate buffers
  ErrorCode  rval = mb-> tag_get_bytes(tag_handle,  bytes_per_tag) ;
  if (rval!=MB_SUCCESS) return rval;

  bool specified_ids = send_IDs_map.size() > 0;
  int indexReq=0;
  if (!specified_ids) // original send
  {
    std::map<int, Range> split_ranges;
    rval = split_owned_range_general ( owned, split_ranges);
    if (rval!=MB_SUCCESS) return rval;

    // use the buffers data structure to allocate memory for sending the tags
    sendReqs.resize(split_ranges.size());

    for (std::map<int, Range>::iterator it=split_ranges.begin(); it!=split_ranges.end(); it++)
    {
      int receiver_proc = it->first;
      Range ents = it->second; // primary entities, with the tag data
      int size_buffer = 4 + bytes_per_tag*(int)ents.size(); // hopefully, below 2B; if more, we have a big problem ...
      ParallelComm::Buffer * buffer = new ParallelComm::Buffer(size_buffer);

      buffer->reset_ptr(sizeof(int));
      // copy tag data to buffer->buff_ptr, and send the buffer (we could have used regular char arrays)
      rval = mb->tag_get_data(tag_handle, ents, (void*)(buffer->buff_ptr) );
      if (rval!=MB_SUCCESS) return rval;
      *((int*)buffer->mem_ptr) = size_buffer;
      //int size_pack = buffer->get_current_size(); // debug check
      ierr = MPI_Isend(buffer->mem_ptr, size_buffer, MPI_CHAR, receiver_proc, 222, jcomm, &sendReqs[indexReq]); // we have to use global communicator
      if (ierr!=0) return MB_FAILURE;
      indexReq++;
      localSendBuffs.push_back(buffer); // we will release them after nonblocking sends are completed
    }
  }
  else
  {
    // we know that we will need to send some tag data in a specific order
    // first, get the ids of the local elements, from owned Range; arrange the buffer in order of increasing global id
    Tag gidTag;
    rval = mb->tag_get_handle("GLOBAL_ID", gidTag); MB_CHK_ERR ( rval );
    std::vector<int> gids;
    gids.resize(owned.size());
    rval = mb->tag_get_data(gidTag, owned, &gids[0]);  MB_CHK_ERR ( rval );
    std::map<int, EntityHandle> gidToHandle;
    size_t i=0;
    for (Range::iterator it=owned.begin(); it!= owned.end(); it++)
    {
      EntityHandle eh = *it;
      gidToHandle[gids[i++]] = eh;
    }
    // now, pack the data and send it
    for (std::map<int ,std::vector<int> >::iterator mit =send_IDs_map.begin(); mit!=send_IDs_map.end(); mit++)
    {
      int receiver_proc = mit->first;
      std::vector<int> & eids = mit->second;
      int size_buffer = 4 + bytes_per_tag*(int)eids.size(); // hopefully, below 2B; if more, we have a big problem ...
      ParallelComm::Buffer * buffer = new ParallelComm::Buffer(size_buffer);
      buffer->reset_ptr(sizeof(int));
      // copy tag data to buffer->buff_ptr, and send the buffer (we could have used regular char arrays)
      for (std::vector<int>::iterator it=eids.begin(); it!=eids.end(); it++)
      {
        int eID = *it;
        EntityHandle eh = gidToHandle[eID];

        rval = mb->tag_get_data(tag_handle, &eh, 1, (void*)(buffer->buff_ptr) );  MB_CHK_ERR ( rval );
        buffer->buff_ptr+=bytes_per_tag;
      }
      *((int*)buffer->mem_ptr) = size_buffer;
        //int size_pack = buffer->get_current_size(); // debug check
      ierr = MPI_Isend(buffer->mem_ptr, size_buffer, MPI_CHAR, receiver_proc, 222, jcomm, &sendReqs[indexReq]); // we have to use global communicator
      if (ierr!=0) return MB_FAILURE;
      indexReq++;
      localSendBuffs.push_back(buffer); // we will release them after nonblocking sends are completed
    }
  }

  return MB_SUCCESS;
}

ErrorCode ParCommGraph::receive_tag_values (MPI_Comm jcomm, ParallelComm *pco, Range & owned, Tag & tag_handle )
{
  // opposite to sending, we will use blocking receives
  int ierr;
  MPI_Status status;
  // basically, owned.size() needs to be equal to sum(corr_sizes)
  // get info about the tag size, type, etc
  Core * mb = (Core*)pco->get_moab();
  // get info about the tag
  //! Get the size of the specified tag in bytes
  int bytes_per_tag; // we need to know, to allocate buffers
  ErrorCode  rval = mb-> tag_get_bytes(tag_handle,  bytes_per_tag) ;
  if (rval!=MB_SUCCESS) return rval;

  bool specified_ids = recv_IDs_map.size() > 0;
  if (!specified_ids)
  {
    std::map<int, Range> split_ranges;
    rval = split_owned_range_general ( owned, split_ranges);
    if (rval!=MB_SUCCESS) return rval;

    // use the buffers data structure to allocate memory for sending the tags
    for (std::map<int, Range>::iterator it=split_ranges.begin(); it!=split_ranges.end(); it++)
    {
      int sender_proc = it->first;
      Range ents = it->second; // primary entities, with the tag data, we will receive
      int size_buffer = 4 + bytes_per_tag*(int)ents.size(); // hopefully, below 2B; if more, we have a big problem ...
      ParallelComm::Buffer * buffer = new ParallelComm::Buffer(size_buffer);

      buffer->reset_ptr(sizeof(int));

      *((int*)buffer->mem_ptr) = size_buffer;
      //int size_pack = buffer->get_current_size(); // debug check

      ierr = MPI_Recv (buffer->mem_ptr, size_buffer, MPI_CHAR, sender_proc, 222, jcomm, &status);
      if (ierr!=0) return MB_FAILURE;
      // now set the tag
      // copy to tag
      rval = mb->tag_set_data(tag_handle, ents, (void*)(buffer->mem_ptr + 4) );
      delete buffer; // no need for it afterwards
      if (rval!=MB_SUCCESS) return rval;

    }
  }
  else // receive buffer, then extract tag data
  {
    // we know that we will need to receive some tag data in a specific order (by ids stored)
    // first, get the ids of the local elements, from owned Range; unpack the buffer in order
    Tag gidTag;
    rval = mb->tag_get_handle("GLOBAL_ID", gidTag); MB_CHK_ERR ( rval );
    std::vector<int> gids;
    gids.resize(owned.size());
    rval = mb->tag_get_data(gidTag, owned, &gids[0]);  MB_CHK_ERR ( rval );
    std::map<int, EntityHandle> gidToHandle;
    size_t i=0;
    for (Range::iterator it=owned.begin(); it!= owned.end(); it++)
    {
      EntityHandle eh = *it;
      gidToHandle[gids[i++]] = eh;
    }
    //
    // now, unpack the data and set it to the tag
    for (std::map<int ,std::vector<int> >::iterator mit =recv_IDs_map.begin(); mit!=recv_IDs_map.end(); mit++)
    {
      int sender_proc = mit->first;
      std::vector<int> & eids = mit->second;
      int size_buffer = 4 + bytes_per_tag*(int)eids.size(); // hopefully, below 2B; if more, we have a big problem ...
      ParallelComm::Buffer * buffer = new ParallelComm::Buffer(size_buffer);
      buffer->reset_ptr(sizeof(int));
      *((int*)buffer->mem_ptr) = size_buffer;

      // receive the buffer
      ierr = MPI_Recv (buffer->mem_ptr, size_buffer, MPI_CHAR, sender_proc, 222, jcomm, &status);
      if (ierr!=0) return MB_FAILURE;

      // copy tag data to buffer->buff_ptr, and send the buffer (we could have used regular char arrays)
      for (std::vector<int>::iterator it=eids.begin(); it!=eids.end(); it++)
      {
        int eID = *it;
        EntityHandle eh = gidToHandle[eID];

        rval = mb->tag_set_data(tag_handle, &eh, 1, (void*)(buffer->buff_ptr) );  MB_CHK_ERR ( rval );
        buffer->buff_ptr+=bytes_per_tag;
      }
    }

  }
  return MB_SUCCESS;
}

ErrorCode ParCommGraph::distribute_sender_graph_info(MPI_Comm senderComm, std::vector<int> &sendingInfo )
{
  // only the root of the sender has all info
  // it will use direct send/receives, for number of elems to be sent to each receiver, in an array
  if (rootSender)
  {
    // just append the local info for the array , and send some data for each receiver
    int nrecv0 = (int) recv_graph[senderTasks[0]].size();
    for (int k=0; k<nrecv0; k++)
    {
      sendingInfo.push_back(recv_graph[senderTasks[0]][k]);
      sendingInfo.push_back(recv_sizes[senderTasks[0]][k]);
    }
    // the rest of the info will be sent for each sender task
    for (int j=1; j< (int) senderTasks.size(); j++)
    {
      std::vector<int> array_to_send;
      // in the sender comm and sender group , rank to send to is j
      for (int k=0; k<(int) recv_graph[senderTasks[j]].size(); k++)
      {
        array_to_send.push_back( recv_graph[senderTasks[j]][k]);
        array_to_send.push_back( recv_sizes[senderTasks[j]][k]);
      }

      int ierr = MPI_Send(&array_to_send[0], (int)array_to_send.size(), MPI_INT, j, 11, senderComm);
      if (MPI_SUCCESS != ierr)
        return MB_FAILURE;
    }
  }
  else
  {
    // expect to receive some info about where to send local data
    // local array it is max the size of 2 times nreceivers;
    int sizeBuffMax = 2 * (int) receiverTasks.size();
    std::vector<int> data(sizeBuffMax);
    MPI_Status status;
    int ierr = MPI_Recv(&data[0], sizeBuffMax, MPI_INT, 0, 11, senderComm, &status);
    if (MPI_SUCCESS != ierr)
      return MB_FAILURE;
    int count;
    // find out how much data we actually got
    MPI_Get_count(&status, MPI_INT, &count);
    for (int k=0; k<count; k++)
      sendingInfo.push_back(data[k]);

  }
  return MB_SUCCESS;
}

ErrorCode ParCommGraph::settle_send_graph(TupleList & TLcovIDs)
{
  // fill send_IDs_map with data
  // will have "receiving proc" and global id of element
  int n = TLcovIDs.get_n();
  for (int i=0; i<n; i++)
  {
    int to_proc= TLcovIDs.vi_wr[2 * i];
    int globalIdElem= TLcovIDs.vi_wr[2 * i + 1 ];
    send_IDs_map[to_proc].push_back(globalIdElem);
  }
#ifdef VERBOSE
  for (std::map<int ,std::vector<int> >::iterator mit=send_IDs_map.begin(); mit!=send_IDs_map.end(); mit++)
  {
    std::cout <<" towards task " << mit->first << " send: " << mit->second.size() << " cells "<< std::endl;
    for (size_t i=0; i< mit->second.size() ; i++)
    {
      std::cout << " " <<  mit->second[i] ;
    }
    std::cout<<std::endl;
  }
#endif
  return MB_SUCCESS;
}

// this will set recv_IDs_map will store all ids to be received from one sender task
void ParCommGraph::SetReceivingAfterCoverage(std::map<int, std::set<int> > & idsFromProcs) // will make sense only on receivers, right now after cov
{
  for (std::map<int, std::set<int> >::iterator mt=idsFromProcs.begin(); mt!=idsFromProcs.end(); mt++)
  {
    int fromProc = mt->first;
    std::set<int> & setIds = mt->second;
    recv_IDs_map[fromProc].resize( setIds.size() );
    std::vector<int>  & listIDs = recv_IDs_map[fromProc];
    size_t indx = 0;
    for ( std::set<int>::iterator st= setIds.begin(); st!=setIds.end(); st++)
    {
      int valueID = *st;
      listIDs[indx++]=valueID;
    }
  }
  return;
}
} // namespace moab
