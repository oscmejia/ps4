#include "Mesh.hpp"
#include "CS207/Util.hpp"



typedef Mesh<int> MeshType;
typedef MeshType::Node node_type;

using namespace std;
using namespace CS207;

static unsigned fail_count = 0;

template <typename T>
void sf_print(T a, std::string msg = "") {
  (void) a;
  std::cerr << msg << " [Success]" << std::endl;
}

void sf_print(bool sf, std::string msg = "") {
  if (sf)
    std::cerr << msg << " [Success]" << std::endl;
  else {
    std::cerr << msg << " [FAIL]" << std::endl;
    ++fail_count;
  }
}


int main(int argc, char** argv)
{
  (void) argc;
  (void) argv;

  Point p1(CS207::random(), CS207::random(), CS207::random());
  Point p2(CS207::random(), CS207::random(), CS207::random());
  Point p3(CS207::random(), CS207::random(), CS207::random());


  std::vector<node_type> nodes;
  // Define an empty Mesh
  MeshType m;
  nodes.push_back(m.add_node(p1));
  sf_print(m.num_nodes() == 1, "Mesh has 1 Node");
  
  nodes.push_back(m.add_node(p2));
  sf_print(m.num_nodes() == 2, "Mesh has 2 Nodes");
  
  nodes.push_back(m.add_node(p3));
  sf_print(m.num_nodes() == 3, "Mesh has 3 Nodes");
  
  sf_print(m.num_triangles() == 0, "Mesh has no triangles");
  m.add_triangle( nodes[0], nodes[1], nodes[2] );
  sf_print(m.num_triangles() == 1, "Mesh has 1 triangle");





  if (fail_count) {
    std::cerr << "\n" << fail_count
        << (fail_count > 1 ? " FAILURES" : " FAILURE") << std::endl;
    return 1;
  } 
  else
    return 0;

}
