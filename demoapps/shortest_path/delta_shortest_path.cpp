
#include <graphlab.hpp>


#include <stdio.h>
#include <unistd.h>



#include <graphlab/macros_def.hpp>



struct vertex_data {
  float dist;
  vertex_data() : dist(std::numeric_limits<float>::max()) { }
}; // end of vertex_data

struct edge_data { 
  float dist; 
  edge_data() : dist(1) { } 
}; // end of edge data

typedef graphlab::graph<vertex_data, edge_data> graph_type;


class shortest_path_update : public
graphlab::iupdate_functor<graph_type, shortest_path_update> {

private:
  float dist;

public:

  shortest_path_update() : dist(0) { }
  shortest_path_update(const float& dist) : dist(dist) { }
  
  void operator+=(const shortest_path_update& other){
    dist = std::min(dist, other.dist);
  }
  
  consistency_model consistency() const {
    return graphlab::VERTEX_CONSISTENCY;
  }
  
  void operator()(icontext_type& context) {
  
    vertex_data& vdata = context.vertex_data();
    
    if (dist < vdata.dist){  
      vdata.dist = dist;
      // reschedule any affected neighbors
      foreach(edge_type edge, context.out_edges()) {
        const vertex_id_type nbr_id = edge.target();
        const edge_data& edata = context.const_edge_data(edge);
        context.schedule(nbr_id, shortest_path_update(edata.dist + dist));    
      }
    }
    
  }
  
}; // end of shortest path update functor

/**
 * This aggregator sums the distances and finds the longest distance.
 */       
class aggregator :
  public graphlab::iaggregator<graph_type, shortest_path_update, aggregator> {
private:
  float max_dist;
  float sum_dist;
  vertex_id_type furthest_vertex;
public:
  aggregator() : max_dist(0), sum_dist(0), furthest_vertex(0) { }
  void operator()(icontext_type& context) {
    float dist = context.const_vertex_data().dist;
    sum_dist += dist;
    if (dist > max_dist){
      max_dist = dist;
      furthest_vertex = context.vertex_id();
    }
  } // end of operator()
  void operator+=(const aggregator& other) {
    sum_dist += other.sum_dist;
    if (other.max_dist > max_dist){
      max_dist = other.max_dist;
      furthest_vertex = other.furthest_vertex;
    }
  }
  void finalize(iglobal_context_type& context) {
    std::cout << "Total Distance:\t\t" << sum_dist << std::endl 
              << "Longest Distance:\t" <<  max_dist << std::endl
              << "Furthest Vertex:\t"  << furthest_vertex << std::endl;
  }
}; // end of aggregator


int main(int argc, char** argv) {
  // Parse input
  graphlab::command_line_options clopts("Run Shortest Path Algorithm.");
  std::string graph_file;
  std::string format = "snap";
  size_t root = 0;
  clopts.attach_option("graph", &graph_file, graph_file,
                       "The graph file.");
  clopts.attach_option("format", &format, format,
                       "File format.");
  clopts.attach_option("root", &root, root,
                       "The root vertex.");
  clopts.add_positional("graph");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing input." << std::endl;
    return EXIT_FAILURE;
  }

  // Create a graphlab core
  graphlab::core<graph_type, shortest_path_update> core;
  core.set_options(clopts);
  
  std::cout << "Loading graph from file" << std::endl;
  const bool success = graphlab::graph_ops::
    load_structure (graph_file, format, core.graph());

  if(!success) {
    std::cout << "Error in reading file: " << graph_file
              << std::endl;
    return EXIT_FAILURE;
  }

  // Set the root
  core.schedule(root, shortest_path_update(0));
  std::cout << "Running!" << std::endl;
  const double runtime = core.start();
  std::cout << "Runtime:      " << runtime << std::endl;
  std::cout << "Update Count: " << core.last_update_count() << std::endl;
  std::cout << "Update Freq:  " << (core.last_update_count() / runtime) 
            << std::endl;
            
  core.add_aggregator("aggregator", aggregator(), 0);
  core.aggregate_now("aggregator");

}







