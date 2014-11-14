#include "Mesh.hpp"
#include "CS207/Util.hpp"



typedef Mesh<int,int,int> MeshType;
typedef MeshType::Node node_type;
typedef MeshType::Triangle triangle_type;

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
  Point p4(CS207::random(), CS207::random(), CS207::random());
  Point p5(CS207::random(), CS207::random(), CS207::random());
  Point p6(CS207::random(), CS207::random(), CS207::random());


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
  auto _t1 = m.add_triangle( nodes[0], nodes[1], nodes[2] );
  sf_print(m.num_triangles() == 1, "Mesh has 1 triangle");



  auto first = m.triangle_begin();
  triangle_type t1 = *first;

  sf_print(t1.index() == 0, "Triangle has indx 0");

  sf_print(t1.node(0).position() == p1 || t1.node(0).position() == p2 || t1.node(0).position() == p3, "p1 is stored in one of the nodes in triangle1" );
  sf_print(t1.node(1).position() == p1 || t1.node(1).position() == p2 || t1.node(1).position() == p3, "p2 is stored in one of the nodes in triangle1"  );
  sf_print(t1.node(2).position() == p1 || t1.node(2).position() == p2 || t1.node(2).position() == p3, "p3 is stored in one of the nodes in triangle1"  );
  sf_print(t1.node(2).position() != p4 && t1.node(2).position() != p4 && t1.node(2).position() != p4, "p4 should not be one of the nodes in triangle1"  );
  
  sf_print(m.triangle(0) == _t1, "get triangle at position 0 and compare with _t1");

  nodes.push_back(m.add_node(p4));
  auto _t2 = m.add_triangle( nodes[0], nodes[1], nodes[3] );
  sf_print(m.num_triangles() == 2, "Mesh has 2 triangle");
  sf_print(m.num_nodes() == 4, "Mesh has 4 Nodes");
  sf_print(m.triangle(0) != _t2, "get triangle at position 0 and compare with _t2");
  sf_print(m.triangle(1) == _t2, "get triangle at position 1 and compare with _t2");

  nodes.push_back(m.add_node(p5));
  nodes.push_back(m.add_node(p6));
  sf_print(m.num_triangles() == 2, "Mesh has 2 triangle");
  sf_print(m.num_nodes() == 6, "Mesh has 6 Nodes");

  auto _t3 = m.add_triangle( nodes[1], nodes[3], nodes[4] );
  sf_print(m.num_triangles() == 3, "Mesh has 3 triangle");

  auto _t4 = m.add_triangle( nodes[3], nodes[4], nodes[5] );
  sf_print(m.num_triangles() == 4, "Mesh has 4 triangle");
  sf_print(m.triangle(2) == _t3, "get triangle at position 2 and compare with _t3");
  sf_print(m.triangle(3) == _t4, "get triangle at position 3 and compare with _t4");
  
  ++first;
  triangle_type t2 = *first;
  sf_print(t2.node(0).position() == p1 || t2.node(0).position() == p2 || t2.node(0).position() == p4, "p1 is stored in one of the nodes in triangle2" );
  sf_print(t2.node(1).position() == p1 || t2.node(1).position() == p2 || t2.node(1).position() == p4, "p2 is stored in one of the nodes in triangle2"  );
  sf_print(t2.node(2).position() == p1 || t2.node(2).position() == p2 || t2.node(2).position() == p4, "p4 is stored in one of the nodes in triangle2"  );
  sf_print(t2.node(2).position() != p3 && t2.node(2).position() != p3 && t2.node(2).position() != p3, "p3 should not be one of the nodes in triangle2"  );
  sf_print(t2.node(2).position() != p5 && t2.node(2).position() != p5 && t2.node(2).position() != p5, "p5 should not be one of the nodes in triangle2"  );
  sf_print(t2.node(2).position() != p6 && t2.node(2).position() != p6 && t2.node(2).position() != p6, "p6 should not be one of the nodes in triangle2"  );

  ++first;
  triangle_type t3 = *first;
  sf_print(t3.node(0).position() == p5 || t3.node(0).position() == p2 || t3.node(0).position() == p4, "p2 is stored in one of the nodes in triangle3" );
  sf_print(t3.node(1).position() == p5 || t3.node(1).position() == p2 || t3.node(1).position() == p4, "p4 is stored in one of the nodes in triangle3"  );
  sf_print(t3.node(2).position() == p5 || t3.node(2).position() == p2 || t3.node(2).position() == p4, "p5 is stored in one of the nodes in triangle3"  );
  sf_print(t3.node(2).position() != p1 && t3.node(2).position() != p1 && t3.node(2).position() != p1, "p1 should not be one of the nodes in triangle3"  );
  sf_print(t3.node(2).position() != p3 && t3.node(2).position() != p3 && t3.node(2).position() != p3, "p3 should not be one of the nodes in triangle3"  );
  sf_print(t3.node(2).position() != p6 && t3.node(2).position() != p6 && t3.node(2).position() != p6, "p6 should not be one of the nodes in triangle3"  );
  
  ++first;
  triangle_type t4 = *first;
  sf_print(t4.node(0).position() == p4 || t4.node(0).position() == p5 || t4.node(0).position() == p6, "p4 is stored in one of the nodes in triangle4" );
  sf_print(t4.node(1).position() == p4 || t4.node(1).position() == p5 || t4.node(1).position() == p6, "p5 is stored in one of the nodes in triangle4"  );
  sf_print(t4.node(2).position() == p4 || t4.node(2).position() == p5 || t4.node(2).position() == p6, "p6 is stored in one of the nodes in triangle4"  );
  sf_print(t4.node(2).position() != p1 && t4.node(2).position() != p1 && t4.node(2).position() != p1, "p1 should not be one of the nodes in triangle4"  );
  sf_print(t4.node(2).position() != p2 && t4.node(2).position() != p2 && t4.node(2).position() != p2, "p2 should not be one of the nodes in triangle4"  );
  sf_print(t4.node(2).position() != p3 && t4.node(2).position() != p3 && t4.node(2).position() != p3, "p3 should not be one of the nodes in triangle4"  );

  sf_print(t4.index() == 3, "triangle4 has idx 3" );
  sf_print(t3.index() == 2, "triangle3 has idx 2" );
  sf_print(t2.index() == 1, "triangle2 has idx 1" );
  sf_print(t1.index() == 0, "triangle1 has idx 0" );

  sf_print(t1 != t2, "triangles are different" );
  sf_print(t1 != t3, "triangles are different" );
  sf_print(t1 != t4, "triangles are different" );
  sf_print(t2 != t3, "triangles are different" );
  sf_print(t2 != t4, "triangles are different" );
  sf_print(t3 != t4, "triangles are different" );

  auto other_first = m.triangle_begin();
  triangle_type other_t1 = *other_first;
  sf_print(t1 == other_t1, "triangles are equal" );

  auto f = m.triangle_begin();
  auto l = m.triangle_end();
  int tri_found = 0;
  while(f != l){
    ++tri_found;
    ++f;
  }
  sf_print(tri_found == 4, "triangle iterator found 4 triangles" );



  
  
  



  if (fail_count) {
    std::cerr << "\n" << fail_count
        << (fail_count > 1 ? " FAILURES" : " FAILURE") << std::endl;
    return 1;
  } 
  else
    return 0;

}
