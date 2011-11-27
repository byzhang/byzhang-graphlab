/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


/**
 * Functionality: The code solves the linear system Ax = b using
 * The Jacobi algorithm. (A is a square matrix). 
 * A assumed to be full column rank.  Algorithm is described
 * http://en.wikipedia.org/wiki/Jacobi_method
 * Written by Danny Bickson 
 */

#ifndef JACOBI_HPP
#define JACOBI_HPP

#include <cmath>
#include <cstdio>
#include <limits>
#include <iostream>
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
using namespace graphlab;

#include <graphlab/macros_def.hpp>

bool debug = false;
int max_iter = 50;
int * stats_per_round = NULL;
int iiter = 0; //current iteration

enum kcore_output_fields{
  KCORE_INDEX = 1
};

struct vertex_data {
  bool active;
  int kcore, degree;

  vertex_data() : active(true), kcore(-1), degree(0)  {}

  void add_self_edge(double value) { }

  void set_val(double value, int field_type) { 
  }
  //only one output for jacobi - solution x
  double get_output(int field_type){ 
    if (field_type == KCORE_INDEX)
      return kcore;
    return -1;
  }
}; // end of vertex_data

struct edge_data {
  edge_data()  { }
  //compatible with parser which have edge value (we don't need it)
  edge_data(double val)  { }
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;

void calc_initial_degree(graph_type * g, matrix_descriptor & desc){
  for (int i=0; i< desc.total(); i++){
     vertex_data & data = g->vertex_data(i);
     data.degree = g->out_edge_ids(i).size() + g->in_edge_ids(i).size();
  }
}


struct kcore_update :
  public graphlab::iupdate_functor<graph_type, kcore_update> {
  void operator()(icontext_type& context) {
  } 
};

class accumulator :
  public graphlab::iaccumulator<graph_type, kcore_update, accumulator> {
private:
  int num_active;
public:
  accumulator() : num_active(0) { }

  void operator()(icontext_type& context) {
    vertex_data & vdata = context.vertex_data();
    int cur_iter = iiter + 1;
    if (vdata.degree <= cur_iter){
       vdata.active = false;
       vdata.kcore = cur_iter;
       vdata.degree = 0;
    } 
    int links= 0;
    edge_list_type outedgeid = context.out_edge_ids();
    edge_list_type inedgeid = context.in_edge_ids();

    for(size_t i = 0; i < outedgeid.size(); i++) {
        const vertex_data & other = context.const_vertex_data(context.target(outedgeid[i]));
        if (other.active)
	  links++;
    }
    for (size_t i =0; i < inedgeid.size(); i++){
	const vertex_data & other = context.const_vertex_data(context.source(inedgeid[i]));
        if (other.active)
          links++;
    }
    if (links <= cur_iter){
        vdata.active = false;
        vdata.kcore = cur_iter;
        vdata.degree = 0;
    }
    else {
        vdata.degree = links;
    }
    if (vdata.active)
      num_active++;
  };

  void operator+=(const accumulator& other) { 
    num_active += other.num_active;
  }

  void finalize(iglobal_context_type& context) {
   printf("Number of active nodes in round %d is %d\n", iiter, num_active);
   stats_per_round[iiter] = num_active;
   if (iiter >= 2)
    printf("Nodes removed from core in this round %d\n", stats_per_round[iiter] - stats_per_round[iiter-1]);

   if (num_active == 0)
     context.terminate(); 
   }
}; // end of  accumulator





int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string datafile;
  std::string format = "matrixmarket";
  int unittest = 0;

  clopts.attach_option("data", &datafile, datafile,
                       "matrix A input file");
  clopts.add_positional("data");
  clopts.attach_option("format", &format, format, "matrix format");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=TBD");
  clopts.attach_option("max_iter", &max_iter, max_iter, "maximal number of cores");

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  stats_per_round = new int[max_iter];

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab Linear solver library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl 
    << "Currently implemented algorithms are: Gaussian Belief Propagation, "
    << "Jacobi method, Conjugate Gradient" << std::endl;



  // Create a core
  graphlab::core<graph_type, kcore_update> core;
  core.set_options(clopts); // Set the engine options

  //unit testing
  if (unittest == 1){
  }

  std::cout << "Load graph" << std::endl;
  matrix_descriptor matrix_info;
  load_graph(datafile, format, matrix_info, core.graph());

  //std::cout << "Schedule all vertices" << std::endl;
  //core.schedule_all(jacobi_update());
 
  accumulator acum;
  core.add_sync("sync", acum, 1000);
  core.add_global("NUM_ACTIVE", int(0));

  graphlab::timer mytimer; mytimer.start();

  for (iiter=0; iiter< max_iter; iiter++){
    core.sync_now("sync");
  }
 
  std::cout << "Jacobi finished in " << mytimer.current_time() << std::endl;

  vec ret = fill_output(&core.graph(), matrix_info, KCORE_INDEX);
  write_output_vector(datafile + "x.out", format, ret);


  if (unittest == 1){
  }

   return EXIT_SUCCESS;
}



#include <graphlab/macros_undef.hpp>
#endif
