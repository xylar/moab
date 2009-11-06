/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

// MOAB performance tests building mapped mesh with nodes and
// hexes created one at a time.  This also creates the node to hex adjacencies.
 
// Different platforms follow different conventions for usage
#if !defined(_MSC_VER) && !defined(__MINGW32__)
#include <sys/resource.h>
#endif
#ifdef SOLARIS
extern "C" int getrusage(int, struct rusage *);
#ifndef RUSAGE_SELF
#include </usr/ucbinclude/sys/rusage.h>
#endif
#endif

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>
#include "MBCore.hpp"
#include "MBReadUtilIface.hpp"
#include "ScdVertexData.hpp"
#include "VertexSequence.hpp"
#include "StructuredElementSeq.hpp"
#include "EntitySequence.hpp"
#include "SequenceManager.hpp"
#include "HomXform.hpp"
#include "MBSkinner.hpp"

double LENGTH = 1.0;
const int DEFAULT_INTERVALS = 50;


void testA(const int nelem);
void testB(const int nelem);
void testC(const int nelem);
void skin_time(const int intervals);
void skin_adj(const int intervals);
void skin_vert(const int intervals);
void print_time(const bool print_em, double &tot_time, double &utime, double &stime,
                double &mem);
void query_vert_to_elem(MBInterface*);
void query_elem_to_vert(MBInterface*);
void query_struct_elem_to_vert(MBInterface*);
void get_time_mem(double &tot_time, double &user_time,
                  double &sys_time, double &tot_mem);
void bulk_construct_mesh( MBInterface* gMB, const int nelem );

void compute_edge(double *start, const int nelem,  const double xint,
                  const int stride) 
{
  for (int i = 1; i < nelem; i++) {
    start[i*stride] = start[0]+i*xint;
    start[nelem+1+i*stride] = start[nelem+1]+i*xint;
    start[2*(nelem+1)+i*stride] = start[2*(nelem+1)]+i*xint;
  }
}

void compute_face(double *a, const int nelem,  const double xint,
                  const int stride1, const int stride2) 
{
    // 2D TFI on a face starting at a, with strides stride1 in ada and stride2 in tse
  for (int j = 1; j < nelem; j++) {
    double tse = j * xint;
    for (int i = 1; i < nelem; i++) {
      double ada = i * xint;
      
      a[i*stride1+j*stride2] = (1.0 - ada)*a[i*stride1]
        + ada*a[i*stride1+nelem*stride2]
        + (1.0 - tse)*a[j*stride2]
        + tse*a[j*stride2+nelem*stride1]
        - (1.0 - tse)*(1.0 - ada)*a[0]
        - (1.0 - tse)*ada*a[nelem*stride1]
        - tse*(1.0 - ada)*a[nelem*stride2]
        - tse*ada*a[nelem*(stride1+stride2)];
      a[nelem+1+i*stride1+j*stride2] = (1.0 - ada)*a[nelem+1+i*stride1]
        + ada*a[nelem+1+i*stride1+nelem*stride2]
        + (1.0 - tse)*a[nelem+1+j*stride2]
        + tse*a[nelem+1+j*stride2+nelem*stride1]
        - (1.0 - tse)*(1.0 - ada)*a[nelem+1+0]
        - (1.0 - tse)*ada*a[nelem+1+nelem*stride1]
        - tse*(1.0 - ada)*a[nelem+1+nelem*stride2]
        - tse*ada*a[nelem+1+nelem*(stride1+stride2)];
      a[2*(nelem+1)+i*stride1+j*stride2] = (1.0 - ada)*a[2*(nelem+1)+i*stride1]
        + ada*a[2*(nelem+1)+i*stride1+nelem*stride2]
        + (1.0 - tse)*a[2*(nelem+1)+j*stride2]
        + tse*a[2*(nelem+1)+j*stride2+nelem*stride1]
        - (1.0 - tse)*(1.0 - ada)*a[2*(nelem+1)+0]
        - (1.0 - tse)*ada*a[2*(nelem+1)+nelem*stride1]
        - tse*(1.0 - ada)*a[2*(nelem+1)+nelem*stride2]
        - tse*ada*a[2*(nelem+1)+nelem*(stride1+stride2)];
    }
  }
}

void build_coords(const int nelem, std::vector<double>& coords) 
{
  double ttime0, ttime1, utime1, stime1, mem1;
  print_time(false, ttime0, utime1, stime1, mem1);
  std::cout << "Before TFI, total memory = " << mem1 / 1.0e6 << " MB." << std::endl;
  
    // allocate the memory
  int numv = nelem+1;
  int numv_sq = numv*numv;
  int tot_numv = numv*numv*numv;
  coords.resize(3*tot_numv);

// use FORTRAN-like indexing
#define VINDEX(i,j,k) (i + (j*numv) + (k*numv_sq))
  int idx;
  double scale1, scale2, scale3;
    // use these to prevent optimization on 1-scale, etc (real map wouldn't have
    // all these equal)
  scale1 = LENGTH/nelem;
  scale2 = LENGTH/nelem;
  scale3 = LENGTH/nelem;

#ifdef REALTFI
    // use a real TFI xform to compute coordinates
    // compute edges
    // i (stride=1)
  compute_edge(&coords[VINDEX(0,0,0)], nelem, scale1, 1);
  compute_edge(&coords[VINDEX(0,nelem,0)], nelem, scale1, 1);
  compute_edge(&coords[VINDEX(0,0,nelem)], nelem, scale1, 1);
  compute_edge(&coords[VINDEX(0,nelem,nelem)], nelem, scale1, 1);
    // j (stride=numv)
  compute_edge(&coords[VINDEX(0,0,0)], nelem, scale1, numv);
  compute_edge(&coords[VINDEX(nelem,0,0)], nelem, scale1, numv);
  compute_edge(&coords[VINDEX(0,0,nelem)], nelem, scale1, numv);
  compute_edge(&coords[VINDEX(nelem,0,nelem)], nelem, scale1, numv);
    // k (stride=numv^2)
  compute_edge(&coords[VINDEX(0,0,0)], nelem, scale1, numv_sq);
  compute_edge(&coords[VINDEX(nelem,0,0)], nelem, scale1, numv_sq);
  compute_edge(&coords[VINDEX(0,nelem,0)], nelem, scale1, numv_sq);
  compute_edge(&coords[VINDEX(nelem,nelem,0)], nelem, scale1, numv_sq);

    // compute faces
    // i=0, nelem
  compute_face(&coords[VINDEX(0,0,0)], nelem, scale1, numv, numv_sq);
  compute_face(&coords[VINDEX(nelem,0,0)], nelem, scale1, numv, numv_sq);
    // j=0, nelem
  compute_face(&coords[VINDEX(0,0,0)], nelem, scale1, 1, numv_sq);
  compute_face(&coords[VINDEX(0,nelem,0)], nelem, scale1, 1, numv_sq);
    // k=0, nelem
  compute_face(&coords[VINDEX(0,0,0)], nelem, scale1, 1, numv);
  compute_face(&coords[VINDEX(0,0,nelem)], nelem, scale1, 1, numv);

    // initialize corner indices
  int i000 = VINDEX(0,0,0);
  int ia00 = VINDEX(nelem,0,0);
  int i0t0 = VINDEX(0,nelem,0);
  int iat0 = VINDEX(nelem,nelem,0);
  int i00g = VINDEX(0,0,nelem);
  int ia0g = VINDEX(nelem,0,nelem);
  int i0tg = VINDEX(0,nelem,nelem);
  int iatg = VINDEX(nelem,nelem,nelem);
  double cX, cY, cZ;
  int adaInts = nelem;
  int tseInts = nelem;
  int gammaInts = nelem;
  
  
  for (int i=1; i < nelem; i++) {
    for (int j=1; j < nelem; j++) {
      for (int k=1; k < nelem; k++) {
        idx = VINDEX(i,j,k);
        double tse = i*scale1;
        double ada = j*scale2;
        double gamma = k*scale3;
        double tm1 = 1.0 - tse;
        double am1 = 1.0 - ada;
        double gm1 = 1.0 - gamma;

        cX = gm1 *   (am1*(tm1*coords[i000] + tse*coords[i0t0])  +
                      ada*(tm1*coords[ia00] + tse*coords[iat0])) +
          gamma * (am1*(tm1*coords[i00g] + tse*coords[i0tg])  +
                   ada*(tm1*coords[ia0g] + tse*coords[iatg]));

        cY = gm1 *   (am1*(tm1*coords[i000] + tse*coords[i0t0])  +
                      ada*(tm1*coords[ia00] + tse*coords[iat0])) +
          gamma * (am1*(tm1*coords[i00g] + tse*coords[i0tg])  +
                   ada*(tm1*coords[ia0g] + tse*coords[iatg]));

        cZ = gm1 *   (am1*(tm1*coords[i000] + tse*coords[i0t0])  +
                      ada*(tm1*coords[ia00] + tse*coords[iat0])) +
          gamma * (am1*(tm1*coords[i00g] + tse*coords[i0tg])  +
                   ada*(tm1*coords[ia0g] + tse*coords[iatg]));

        double *ai0k = &coords[VINDEX(k,0,i)];
        double *aiak = &coords[VINDEX(k,adaInts,i)];
        double *a0jk = &coords[VINDEX(k,j,0)];
        double *atjk = &coords[VINDEX(k,j,tseInts)];
        double *aij0 = &coords[VINDEX(0,j,i)];
        double *aijg = &coords[VINDEX(gammaInts,j,i)];
  
        coords[VINDEX(i,j,k)] = (   am1*ai0k[0] 
                                    + ada*aiak[0] 
                                    + tm1*a0jk[0] 
                                    + tse*atjk[0]
                                    + gm1*aij0[0] 
                                    + gamma*aijg[0] )/2.0 - cX/2.0;

        coords[nelem+1+VINDEX(i,j,k)] = (   am1*ai0k[nelem+1] 
                                            + ada*aiak[nelem+1] 
                                            + tm1*a0jk[nelem+1] 
                                            + tse*atjk[nelem+1]
                                            + gm1*aij0[nelem+1] 
                                            + gamma*aijg[nelem+1] )/2.0 - cY/2.0;

        coords[2*(nelem+1)+VINDEX(i,j,k)] = (   am1*ai0k[2*(nelem+1)] 
                                                + ada*aiak[2*(nelem+1)] 
                                                + tm1*a0jk[2*(nelem+1)] 
                                                + tse*atjk[2*(nelem+1)]
                                                + gm1*aij0[2*(nelem+1)] 
                                                + gamma*aijg[2*(nelem+1)] )/2.0 - cZ/2.0;
      }
    }
  }
  

#else
  for (int i=0; i < numv; i++) {
    for (int j=0; j < numv; j++) {
      for (int k=0; k < numv; k++) {
        idx = VINDEX(i,j,k);
          // blocked coordinate ordering
        coords[idx] = i*scale1;
        coords[tot_numv+idx] = j*scale2;
        coords[2*tot_numv+idx] = k*scale3;
      }
    }
  }
#endif
  print_time(false, ttime1, utime1, stime1, mem1);
  std::cout << "MOAB: TFI time = " << ttime1-ttime0 << " sec, memory = " 
            << mem1/1.0e6 << " MB." << std::endl;
}

void build_connect(const int nelem, const MBEntityHandle vstart, MBEntityHandle *connect) 
{
  MBEntityHandle vijk;
  int numv = nelem + 1;
  int numv_sq = numv*numv;
  int idx = 0;
  for (int i=0; i < nelem; i++) {
    for (int j=0; j < nelem; j++) {
      for (int k=0; k < nelem; k++) {
        vijk = vstart+VINDEX(i,j,k);
        connect[idx++] = vijk;
        connect[idx++] = vijk+1;
        connect[idx++] = vijk+1+numv;
        connect[idx++] = vijk+numv;
        connect[idx++] = vijk+numv*numv;
        connect[idx++] = vijk+1+numv*numv;
        connect[idx++] = vijk+1+numv+numv*numv;
        connect[idx++] = vijk+numv+numv*numv;
        assert(i <= numv*numv*numv);
      }
    }
  }
}

typedef void (*test_func_t)( int );
const struct {
  std::string testName;
  test_func_t testFunc;
  std::string testDesc;
} TestList[] = {
 { "struct",    &testA,     "Conn. and adj. query time for structured mesh" },
 { "bulk",      &testB,     "Conn. and adj. query time for bulk-created mesh" },
 { "indiv",     &testC,     "Conn. and adj. query time for per-entity created mesh" },
 { "skin",      &skin_time, "Test time to get skin quads" },
 { "skin_adj",  &skin_adj,  "Test skin quads using vert-to-elem adjacencies" },
 { "skin_vert", &skin_vert, "Test time to get skin verts" },
};
const int TestListSize = sizeof(TestList)/sizeof(TestList[0]); 

void usage( const char* argv0, bool error = true )
{
  std::ostream& str = error ? std::cerr : std::cout;
  str << "Usage: " << argv0 << " -i <ints_per_side> <test_name> [<test_name2> ...]" << std::endl;
  str << "       " << argv0 << " [-h|-l]" << std::endl;
  if (error)
    return;
  str << "  -i  : specify interverals per side (num hex = ints^3, default: " << DEFAULT_INTERVALS << std::endl;
  str << "  -h  : print this help text." << std::endl;
  str << "  -l  : list available tests"  << std::endl;
}

void list_tests( ) 
{
  unsigned max_test_name = 0, max_test_desc = 0;
  for (int i = 0; i < TestListSize; ++i) {
    if (TestList[i].testName.size() > max_test_name)
      max_test_name = TestList[i].testName.size();
    if (TestList[i].testDesc.size() > max_test_desc)
      max_test_desc = TestList[i].testDesc.size();
  }
  std::cout << std::setw(max_test_name) << "NAME" << "   " 
            << std::setw(max_test_desc) << std::left << "DESCRIPTION" << std::endl
            << std::setfill('-') << std::setw(max_test_name) << "" << "   " 
            << std::setfill('-') << std::setw(max_test_desc) << "" 
            << std::setfill(' ') << std::endl;
  for (int i = 0; i < TestListSize; ++i) 
    std::cout << std::setw(max_test_name) << TestList[i].testName << " : "
              << std::setw(max_test_desc) << std::left << TestList[i].testDesc << std::endl;
}  
  
int main(int argc, char* argv[])
{
  int intervals = DEFAULT_INTERVALS;
  std::vector<test_func_t> test_list;
  bool did_help = false;
  bool expect_ints = false;
  for (int i = 1; i < argc; ++i) {
    if (expect_ints) {
      expect_ints = false;
      char* endptr;
      intervals = strtol( argv[i], &endptr, 0 );
      if (intervals < 1) {
        usage(argv[0]);
        std::cerr << "Invalid interval count: " << intervals << std::endl;
        return 1;
      }
    }
    else if (*argv[i] == '-') { // flag
      for (int j = 1; argv[i][j]; ++j) {
        switch (argv[i][j]) {
          case 'i': expect_ints = true; break;
          case 'h': did_help = true; usage(argv[0],false); break;
          case 'l': did_help = true; list_tests(); break;
          default:
            usage(argv[0]);
            std::cerr << "Invalid option: -" << argv[i][j] << std::endl;
            return 1;
        }
      }
    }
    else {
      int j = -1;
      do {
        ++j;
        if (j >= TestListSize) {
          usage(argv[0]);
          std::cerr << "Invalid test name: " << argv[i] << std::endl;
          return 1;
        }
      } while (TestList[j].testName != argv[i]);
      test_list.push_back( TestList[j].testFunc );
    }
  }
  
    // error if no input
  if (test_list.empty() && !did_help) {
    usage(argv[0]);
    std::cerr << "No tests specified" << std::endl;
    return 1;
  }
  
    // now run the tests
  for (std::vector<test_func_t>::iterator i = test_list.begin(); i != test_list.end(); ++i) {
    test_func_t fptr = *i;
    fptr(intervals);
  }
  
  return 0;
}

void query_elem_to_vert(MBInterface* gMB)
{
  MBRange all_hexes;
  MBErrorCode result = gMB->get_entities_by_type(0, MBHEX, all_hexes);
  const MBEntityHandle *connect;
  int num_connect;
  double dum_coords[24];
  for (MBRange::iterator eit = all_hexes.begin(); eit != all_hexes.end(); eit++) {
    result = gMB->get_connectivity(*eit, connect, num_connect);
    assert(MB_SUCCESS == result);
    result = gMB->get_coords(connect, num_connect, dum_coords);
    assert(MB_SUCCESS == result);

      // compute the centroid
    double centroid[3] = {0.0, 0.0, 0.0};
    for (int j = 0; j < 24;) {
      centroid[0] += dum_coords[j++];
      centroid[1] += dum_coords[j++];
      centroid[2] += dum_coords[j++];
    }
  }
}

void query_vert_to_elem(MBInterface* gMB)
{
  MBRange all_verts;
  std::vector<MBEntityHandle> neighbor_hexes;
  MBErrorCode result = gMB->get_entities_by_type(0, MBVERTEX, all_verts);
  assert(MB_SUCCESS == result);
  for (MBRange::iterator vit = all_verts.begin(); vit != all_verts.end(); vit++) {
    neighbor_hexes.clear();
    result = gMB->get_adjacencies(&(*vit), 1, 3, false, neighbor_hexes);
    assert(MB_SUCCESS == result);
  }
}

void query_struct_elem_to_vert(MBInterface* gMB)
{
    // assumes brick mapped mesh with handles starting at zero
  MBRange all_hexes;
  MBErrorCode result = gMB->get_entities_by_type(0, MBHEX, all_hexes);
  double dum_coords[24];
  std::vector<MBEntityHandle> connect;
  for (MBRange::iterator eit = all_hexes.begin(); eit != all_hexes.end(); eit++) {
    result = gMB->get_connectivity(&(*eit), 1, connect);
    assert(MB_SUCCESS == result);
    result = gMB->get_coords(&connect[0], connect.size(), dum_coords);
    assert(MB_SUCCESS == result);
  }
}

void print_time(const bool print_em, double &tot_time, double &utime, double &stime,
                double &tot_mem) 
{
  get_time_mem(tot_time, utime, stime, tot_mem);
  
  if (print_em) {
    std::cout << "User, system, total time = " << utime << ", " << stime 
              << ", " << tot_time << std::endl;
    std::cout << "Total memory = " << tot_mem / 1.0e6 << " MB." << std::endl;
  }
}

#if defined(_MSC_VER) || defined(__MINGW32__)
void get_time_mem(double &tot_time, double &user_time,
                  double &sys_time, double &tot_mem) 
{
  tot_time = 0;
  user_time = (double)clock() / CLOCKS_PER_SEC;
  sys_time = 0;
  tot_mem = 0;
}
#else
void get_time_mem(double &tot_time, double &user_time,
                  double &sys_time, double &tot_mem) 
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  user_time = (double)r_usage.ru_utime.tv_sec +
    ((double)r_usage.ru_utime.tv_usec/1.e6);
  sys_time = (double)r_usage.ru_stime.tv_sec +
    ((double)r_usage.ru_stime.tv_usec/1.e6);
  tot_time = user_time + sys_time;
  tot_mem = 0;
  if (0 != r_usage.ru_maxrss) {
    tot_mem = r_usage.ru_idrss; 
  }
  else {
      // this machine doesn't return rss - try going to /proc
      // print the file name to open
    char file_str[4096], dum_str[4096];
    int file_ptr = -1, file_len;
    file_ptr = open("/proc/self/stat", O_RDONLY);
    file_len = read(file_ptr, file_str, sizeof(file_str)-1);
    if (file_len == 0) return;
    
    close(file_ptr);
    file_str[file_len] = '\0';
      // read the preceeding fields and the ones we really want...
    int dum_int;
    unsigned int dum_uint, vm_size, rss;
    int num_fields = sscanf(file_str, 
                            "%d " // pid
                            "%s " // comm
                            "%c " // state
                            "%d %d %d %d %d " // ppid, pgrp, session, tty, tpgid
                            "%u %u %u %u %u " // flags, minflt, cminflt, majflt, cmajflt
                            "%d %d %d %d %d %d " // utime, stime, cutime, cstime, counter, priority
                            "%u %u " // timeout, itrealvalue
                            "%d " // starttime
                            "%u %u", // vsize, rss
                            &dum_int, 
                            dum_str, 
                            dum_str, 
                            &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
                            &dum_uint, &dum_uint, &dum_uint, &dum_uint, &dum_uint,
                            &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
                            &dum_uint, &dum_uint, 
                            &dum_int,
                            &vm_size, &rss);
    if (num_fields == 24)
      tot_mem = ((double)vm_size);
  }
}
#endif

void testA( int nelem ) 
{
  MBCore moab;
  MBInterface* gMB = &moab;
  double ttime0, ttime1, ttime2, ttime3, utime, stime, mem;
  
  print_time(false, ttime0, utime, stime, mem);
  std::cout << "Ready to read model into MOAB; memory  = " << mem/1.0e6 << " MB." 
            << std::endl;

    // make a 3d block of vertices
  EntitySequence *dum_seq = NULL;
  ScdVertexData *vseq = NULL;
  StructuredElementSeq *eseq = NULL;
  SequenceManager *seq_mgr = moab.sequence_manager();
  HomCoord vseq_minmax[2] = {HomCoord(0,0,0), 
                             HomCoord(nelem, nelem, nelem)};
  MBEntityHandle vstart, estart;
  
  MBErrorCode result = seq_mgr->create_scd_sequence(vseq_minmax[0], vseq_minmax[1],
                                                    MBVERTEX, 1, vstart, dum_seq);
  if (NULL != dum_seq) vseq = dynamic_cast<ScdVertexData*>(dum_seq->data());
  assert (MB_FAILURE != result && vstart != 0 && dum_seq != NULL && vseq != NULL);
    // now the element sequence
  result = seq_mgr->create_scd_sequence(vseq_minmax[0], vseq_minmax[1], 
                                        MBHEX, 1, estart, dum_seq);
  if (NULL != dum_seq) eseq = dynamic_cast<StructuredElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart != 0 && dum_seq != NULL && eseq != NULL);
  
    // only need to add one vseq to this, unity transform
    // trick: if I know it's going to be unity, just input 3 sets of equivalent points
  result = dynamic_cast<ScdElementData*>(eseq->data())
         ->add_vsequence(vseq, vseq_minmax[0], vseq_minmax[0], vseq_minmax[0], 
                               vseq_minmax[0], vseq_minmax[0], vseq_minmax[0]);
  assert(MB_SUCCESS == result);

  std::vector<double> coords;
  build_coords(nelem, coords);

    // set the coordinates of the vertices
  MBEntityHandle handle;
  int i;
  double dumv[3];
  int num_verts = (nelem + 1)*(nelem + 1)*(nelem + 1);
  for (i = 0, handle = vstart; i < num_verts; i++, handle++) {
    dumv[0] = coords[i]; dumv[1] = coords[num_verts+i]; dumv[2] = coords[2*num_verts+i]; 
    result = gMB->set_coords(&handle, 1, dumv);
    assert(MB_SUCCESS == result);
  }

  print_time(false, ttime1, utime, stime, mem);
  std::cout << "Read model into MOAB; memory  = " << mem/1.0e6 << " MB." << std::endl;

    // query the mesh 2 ways
  query_struct_elem_to_vert(gMB);

  print_time(false, ttime2, utime, stime, mem);
  std::cout << "After E-v query; memory  = " << mem/1.0e6 << " MB." << std::endl;

  query_vert_to_elem(gMB);
  
  print_time(false, ttime3, utime, stime, mem);
  std::cout << "After v-E query; memory  = " << mem/1.0e6 << " MB." << std::endl;

  std::cout << "MOAB scd: nelem, construct, e_to_v query, v_to_e query = " 
            << nelem << ", "
            << ttime1-ttime0 << ", " 
            << ttime2-ttime1 << ", " 
            << ttime3-ttime2 << " seconds" 
            << std::endl;
}

void bulk_construct_mesh( MBInterface* gMB, const int nelem )
{
  int num_verts = (nelem + 1)*(nelem + 1)*(nelem + 1);
  int num_elems = nelem*nelem*nelem;
  MBEntityHandle vstart, estart;

  void *ptr = 0;
  // get the read interface
  gMB->query_interface("MBReadUtilIface", &ptr);
  MBReadUtilIface* readMeshIface = static_cast<MBReadUtilIface*>(ptr);

  // create a sequence to hold the node coordinates
  // get the current number of entities and start at the next slot
  std::vector<double*> coord_arrays;
  MBErrorCode result = readMeshIface->get_node_arrays(3, num_verts, 1, vstart, coord_arrays);
  assert(MB_SUCCESS == result && 1 == vstart &&
         coord_arrays[0] && coord_arrays[1] && coord_arrays[2]);
    // memcpy the coordinate data into place
  std::vector<double> coords;
  build_coords(nelem, coords);
  memcpy(coord_arrays[0], &coords[0], sizeof(double)*num_verts);
  memcpy(coord_arrays[1], &coords[num_verts], sizeof(double)*num_verts);
  memcpy(coord_arrays[2], &coords[2*num_verts], sizeof(double)*num_verts);
  
  MBEntityHandle *conn = 0;
  result = readMeshIface->get_element_array(num_elems, 8, MBHEX, 1, estart, conn);
  assert(MB_SUCCESS == result);
  build_connect(nelem, vstart, conn);
  result = readMeshIface->update_adjacencies(estart, num_elems, 8, conn);
  assert(MB_SUCCESS == result);
}

void testB(const int nelem) 
{
  MBCore moab;
  MBInterface* gMB = &moab;
  double ttime0, ttime1, ttime2, ttime3, utime, stime, mem;
  
  print_time(false, ttime0, utime, stime, mem);
  std::cout << "Ready to read model into MOAB; memory = " << mem/1.0e6 << " MB." << std::endl;

  bulk_construct_mesh( gMB, nelem );

  print_time(false, ttime1, utime, stime, mem);
  std::cout << "Read model into MOAB; memory = " << mem/1.0e6 << " MB." << std::endl;

    // query the mesh 2 ways
  query_elem_to_vert(gMB);

  print_time(false, ttime2, utime, stime, mem);
  std::cout << "After E-v query; memory = " << mem/1.0e6 << " MB." << std::endl;

  query_vert_to_elem(gMB);
  
  print_time(false, ttime3, utime, stime, mem);
  std::cout << "After v-E query; memory = " << mem/1.0e6 << " MB." << std::endl;

  std::cout << "MOAB ucd blocked: nelem, construct, e_to_v query, v_to_e query = " 
            << nelem << ", "
            << ttime1-ttime0 << ", " 
            << ttime2-ttime1 << ", " 
            << ttime3-ttime2 << " seconds" 
            << std::endl;
}

void testC(const int nelem) 
{
  MBCore moab;
  MBInterface* gMB = &moab;
  double ttime0, ttime1, ttime2, ttime3, utime, stime, mem;
  
  print_time(false, ttime0, utime, stime, mem);
  std::cout << "Ready to read data into MOAB; memory = " << mem/1.0e6 << " MB." << std::endl;

    // create the vertices; assume we don't need to keep a list of vertex handles, since they'll
    // be created in sequence
  std::vector<double> coords;
  build_coords(nelem, coords);
  int numv = nelem + 1;
  int numv_sq = numv*numv;
  int num_verts = numv*numv*numv;
  double dum_coords[3] = {coords[0], coords[num_verts], coords[2*num_verts]};
  MBEntityHandle vstart;

  MBErrorCode result = gMB->create_vertex(dum_coords, vstart);
  assert(MB_SUCCESS == result && 1 == vstart);

  MBEntityHandle dum_vert, vijk;
  int i;
  for (i = 1; i < num_verts; i++) {
    dum_coords[0] = coords[i]; 
    dum_coords[1] = coords[num_verts+i]; 
    dum_coords[2] = coords[2*num_verts+i];
    result = gMB->create_vertex(dum_coords, dum_vert);
    assert(MB_SUCCESS == result);
  }

  MBEntityHandle dum_conn[8];
  for (i=0; i < nelem; i++) {
    for (int j=0; j < nelem; j++) {
      for (int k=0; k < nelem; k++) {
        vijk = vstart+VINDEX(i,j,k);
        dum_conn[0] = vijk;
        dum_conn[1] = vijk+1;
        dum_conn[2] = vijk+1+numv;
        dum_conn[3] = vijk+numv;
        dum_conn[4] = vijk+numv*numv;
        dum_conn[5] = vijk+1+numv*numv;
        dum_conn[6] = vijk+1+numv+numv*numv;
        dum_conn[7] = vijk+numv+numv*numv;
        result = gMB->create_element(MBHEX, dum_conn, 8, dum_vert);
        assert(MB_SUCCESS == result);
      }
    }
  }

  print_time(false, ttime1, utime, stime, mem);
  std::cout << "Read data into MOAB; memory  = " << mem/1.0e6 << " MB." << std::endl;

    // query the mesh 2 ways
  query_elem_to_vert(gMB);

  print_time(false, ttime2, utime, stime, mem);
  std::cout << "After E-v query; memory  = " << mem/1.0e6 << " MB." << std::endl;

  query_vert_to_elem(gMB);
  
  print_time(false, ttime3, utime, stime, mem);
  std::cout << "After v-E query; memory  = " << mem/1.0e6 << " MB." << std::endl;

  std::cout << "MOAB ucd indiv: nelem, construct, e_to_v query, v_to_e query = " 
            << nelem << ", "
            << ttime1-ttime0 << ", " 
            << ttime2-ttime1 << ", " 
            << ttime3-ttime2 << " seconds" 
            << std::endl;
}

void skin_time( const int nelem ) 
{
  MBCore moab;
  MBInterface* gMB = &moab;
  MBErrorCode rval;

  bulk_construct_mesh( gMB, nelem );

  MBRange skin, hexes;
  rval = gMB->get_entities_by_dimension( 0, 3, hexes );
  assert(MB_SUCCESS == rval); assert(!hexes.empty());

  MBSkinner tool(gMB);
  clock_t t = clock();
  rval = tool.find_skin( hexes, 2, skin, false );
  double d = ((double)(clock() - t))/CLOCKS_PER_SEC;
  assert(MB_SUCCESS == rval);
  
  std::cout << "Got 2D skin in " << d << " seconds w/out vertex-to-element adjacencies" << std::endl;
}

void skin_adj( const int nelem ) 
{
  MBCore moab;
  MBInterface* gMB = &moab;
  MBErrorCode rval;

  bulk_construct_mesh( gMB, nelem );

  MBRange skin, verts, hexes;
    // force creation of adjacencies
  rval = gMB->get_entities_by_dimension( 0, 0, verts );
  assert(MB_SUCCESS == rval); assert(!verts.empty());
  rval = gMB->get_adjacencies( verts, 3, false, hexes, MBInterface::UNION );
  assert(MB_SUCCESS == rval); assert(!hexes.empty());

  MBSkinner tool(gMB);
  clock_t t = clock();
  rval = tool.find_skin( hexes, 2, skin, true );
  double d = ((double)(clock() - t))/CLOCKS_PER_SEC;
  assert(MB_SUCCESS == rval);

  std::cout << "Got 2D skin in " << d << " seconds with vertex-to-element adjacencies" << std::endl;
}

void skin_vert( const int nelem ) 
{
  MBCore moab;
  MBInterface* gMB = &moab;
  MBErrorCode rval;

  bulk_construct_mesh( gMB, nelem );

  MBRange skin, hexes;
  rval = gMB->get_entities_by_dimension( 0, 3, hexes );
  assert(MB_SUCCESS == rval); assert(!hexes.empty());

  MBSkinner tool(gMB);
  clock_t t = clock();
  rval = tool.find_skin( hexes, 0, skin, false );
  double d = ((double)(clock() - t))/CLOCKS_PER_SEC;
  assert(MB_SUCCESS == rval);
  
  std::cout << "Got 0D skin in " << d << " seconds w/out vertex-to-element adjacencies" << std::endl;
}
