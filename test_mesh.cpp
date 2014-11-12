#include <fstream>
#include "CS207/SDLViewer.hpp"
#include "Mesh.hpp"
#include "CS207/Util.hpp"


typedef Mesh<int> MeshType;
typedef Graph<int, int> GraphType;
typedef GraphType::node_type Node;
typedef GraphType::edge_type Edge;

using namespace std;
using namespace CS207;

static unsigned fail_count = 0;

template <typename T>
void sf_print(T a, string msg = "") {
  (void) a;
  cerr << msg << " [Success]" << endl;
}

void sf_print(bool sf, string msg = "") {
  if (sf)
    cerr << msg << " [Success]" << endl;
  else {
    cerr << msg << " [FAIL]" << endl;
    ++fail_count;
  }
}


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE EDGES_FILE\n";
    exit(1);
  }

  // Define an empty Mesh
  MeshType m;


  
  /*
  g.echo();

  // Get the edge length, should be the same for each edge
  double h = g.edge(0).length();

  sf_print(g.num_nodes() == 25, "Graph has 25 Nodes");
  sf_print(g.num_edges() == 40, "Graph has 40 Edges");
  std::cout << "Edges: " << g.num_edges() << std::endl << std::endl;

  std::cout << std::endl << "Remove node 6: " << std::endl;
  sf_print(g.remove_node(g.node(6)), "Remove node 6");
  sf_print(g.num_nodes() == 24, "Graph has 24 Nodes");
  sf_print(g.num_edges() == 36, "Graph has 40 Edges");
  std::cout << "Edges: " << g.num_edges() << " Nodes: " << g.num_nodes() << std::endl << std::endl;

  g.echo();


  int i = 0;
  for (auto it = g.edge_begin(); it != g.edge_end(); ++it) {

    unsigned idx1 = (*it).node1().index();
    unsigned idx2 = (*it).node2().index();

    std::cout << i <<  " - E.node1().index() : " << idx1 << std::endl;
    std::cout << i <<  " - E.node2().index() : " << idx2 << std::endl;

    std::cout << std::endl;
    ++i;
  }
  */

  return 0;
}
