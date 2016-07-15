/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 *
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */




#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif

#include <iostream>
#include "moab/Core.hpp"
#include "TestUtil.hpp"

#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif

#include "Internals.hpp"
#include "TestUtil.hpp"

using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

int main()
{
  Interface* iface0 = new Core;
  Interface* iface1 = new Core;

  const char* filename = 0;

#ifdef MESHDIR
  filename = STRINGIFY(MESHDIR) "/tet_brick.vtk";
#else
#error Specify MESHDIR to compile test
#endif

  ErrorCode err;
  err =  iface0->load_file( filename);CHECK_ERR(err);
  err =  iface1->load_file( filename);CHECK_ERR(err);

  Range tets0, verts0, edges0, edges1, tris0, tris1;
  err = iface0->get_entities_by_dimension(0,3,tets0);CHECK_ERR(err);
  err = iface0->get_adjacencies(tets0,1,true,edges0,Interface::UNION);CHECK_ERR(err);
  err = iface1->get_adjacencies(tets0,1,true,edges1,Interface::UNION);CHECK_ERR(err);
  CHECK_EQUAL(edges0,edges1);
  err = iface0->get_adjacencies(tets0,2,true,tris0,Interface::UNION);CHECK_ERR(err);
  err = iface1->get_adjacencies(tets0,2,true,tris1,Interface::UNION);CHECK_ERR(err);
  CHECK_EQUAL(tris0,tris1);


  std::vector<EntityHandle> conn_squence;

  // At that point iface0 and iface1 should have the same entity handler
  // associated with entities, so we can use one or another, no difference.


  int repeat = 0;
  for(;repeat!=3;repeat++) {

    std::vector<EntityHandle> conn_seq;


    Range to_delete;
    // Build range of quad to delete form iface0
    int ii = 0;
    for(Range::iterator qit=tets0.begin();qit!=tets0.end();qit++,ii++) {
      if(ii%3) {
        to_delete.insert(*qit);
        const EntityHandle* conn;
        int number_nodes = 0;
        iface0->get_connectivity(*qit,conn,number_nodes);
        conn_seq.insert(conn_seq.end(),conn,&conn[number_nodes]);
      }
    }
    // Buidl tange of edges to delete from iface1
    // Create gaps in sequence, to be filled later on.
    for(Range::iterator eit=edges1.begin();eit!=edges1.end();eit++,ii++) {
      if(ii%3) {
        to_delete.insert(*eit);
      }
    }
    for(Range::iterator eit=tris1.begin();eit!=tris1.end();eit++,ii++) {
      if(ii%3) {
        to_delete.insert(*eit);
      }
    }

    err = iface0->delete_entities(to_delete);
    err = iface1->delete_entities(to_delete);

    for(int qq = 0;qq!=conn_seq.size()/4;qq++) {
      EntityHandle q0,q1;
      err = iface1->create_element(MBTET,&conn_seq[4*qq],4,q1);CHECK_ERR(err);
      err = iface0->create_element(MBTET,&conn_seq[4*qq],4,q0);CHECK_ERR(err);
      CHECK(q0==q1);
    }

    tets0.clear();
    err = iface0->get_entities_by_dimension(0,3,tets0);CHECK_ERR(err);
    Range tets1;
    err = iface1->get_entities_by_dimension(0,3,tets1);CHECK_ERR(err);
    CHECK_EQUAL(tets0,tets1);

    // Finally check adjacency, this finally check if code is deterministic
    for(Range::iterator tit = tets1.begin();tit!=tets1.end();tit++) {
      std::vector<EntityHandle> adj0,adj1;
      err = iface0->get_adjacencies(&*tit,1,1,true,adj0);CHECK_ERR(err);
      err = iface1->get_adjacencies(&*tit,1,1,true,adj1);CHECK_ERR(err);
      CHECK_EQUAL(adj0,adj1);
    }
    for(Range::iterator tit = tets1.begin();tit!=tets1.end();tit++) {
      std::vector<EntityHandle> adj0,adj1;
      err = iface0->get_adjacencies(&*tit,1,2,true,adj0);CHECK_ERR(err);
      err = iface1->get_adjacencies(&*tit,1,2,true,adj1);CHECK_ERR(err);
      CHECK_EQUAL(adj0,adj1);
    }

  }

  delete iface0;
  delete iface1;

  return 0;
}
