/**
 * \file gpu_graph_test.cpp   Test gpu_graph on CPU and GPU.
 */

#include <graphlab/gpu/gpu_graph.hpp>

#include <graphlab/macros_def.hpp>

int main(int argc, char* argv[]) {
  using namespace graphlab;

  typedef int vertex_data;
  typedef int edge_data;
  typedef graph<vertex_data, edge_data> cpu_graph_type;
  typedef gpu_graph<vertex_data, edge_data> gpu_graph_type;

  cpu_graph_type cpu_g;
  cpu_g.add_vertex(0);
  cpu_g.add_vertex(1);
  cpu_g.add_vertex(2);
  cpu_g.add_vertex(3);
  cpu_g.add_edge(0, 1, 1);
  cpu_g.add_edge(0, 2, 2);
  cpu_g.add_edge(1, 2, 3);
  cpu_g.add_edge(1, 3, 4);

  return EXIT_SUCCESS;
}

#include <graphlab/macros_undef.hpp>
