/**
 * @file shallow_water.cpp
 * Implementation of a shallow water system using Mesh
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D point list (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
 */

#include <fstream>
#include <cmath>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Point.hpp"
#include "Mesh.hpp"

// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;

/** Water column characteristics */
struct QVar {
  double h;   // Height of column
  double hx;  // Height times average x velocity of column
  double hy;  // Height times average y velocity of column

  /** Default constructor.
   *
   * A default water column is 1 unit high with no velocity. */
  QVar()
      : h(1), hx(0), hy(0) {
  }
  /** Construct the given water column. */
  QVar(double h_, double hx_, double hy_)
      : h(h_), hx(hx_), hy(hy_) {
  }
  // More operators?

  QVar operator+(QVar q){
  	return QVar(h + q.h, hx+q.hx, hy + q.hy);
  };
  
  QVar operator-(QVar q){
  	return QVar(h - q.h, hx - q.hx, hy - q.hy);
  };

  bool operator==(QVar q){
    return (h == q.h && hx == q.hx && hy == q.hy);
  };
  
  QVar operator/(double x){
  	return QVar(h/x, hx/x, hy/x);
  };

  QVar operator*(double x){
    return QVar(h*x, hx*x, hy*x);
  };
  
};

// HW4B: Placeholder for Mesh Type!
// Define NodeData, EdgeData, TriData, etc
// or redefine for your particular Mesh


/** Custom structure of data to store with Triangles */
struct TriData {
  QVar q_bar;  //Qbar for a triangle
  QVar q_tmp;
  std::vector<QVar> edge_fluxes; //debugging purposes, contains the outward flux for each edge 
}; 

/** Custom structure of data to store with Nodes */
struct NodeData {
   QVar q; //Q vector for a node, average of adj triangles
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  QVar total_flux() {
    QVar q;
    for(auto it = fluxes.begin(); it != fluxes.end(); ++it)
      q = q + *it;
    return q;
  }
  std::vector<QVar> fluxes ;  //vector of up to 2 element with fluxes (Qvars) for adj triangles
};

typedef Mesh<TriData, NodeData, EdgeData> MeshType;
typedef MeshType::Node node_type;
typedef MeshType::Edge edge_type;
typedef MeshType::Triangle triangle_type;


/** Function object for calculating shallow-water flux.
 *          |n
 *   T_k    |---> n = (nx,ny)   T_m
 *   QBar_k |                   QBar_m
 *          |
 * @param[in] nx, ny Defines the 2D outward normal vector n = (@a nx, @a ny)
 *            from triangle T_k to triangle T_m. The length of n is equal to the
 *            the length of the edge, |n| = |e|.
 * @param[in] dt The time step taken by the simulation. Used to compute the
 *               Lax-Wendroff dissipation term.
 * @param[in] qk The values of the conserved variables on the left of the edge.
 * @param[in] qm The values of the conserved variables on the right of the edge.
 * @return The flux of the conserved values across the edge e
 */
struct EdgeFluxCalculator {
  QVar operator()(double nx, double ny, double dt,
                  const QVar& qk, const QVar& qm) {
    // Normalize the (nx,ny) vector
    double n_length = std::sqrt(nx*nx + ny*ny);
    nx /= n_length;
    ny /= n_length;

    // The velocities normal to the edge
    double wm = (qm.hx*nx + qm.hy*ny) / qm.h;
    double wk = (qk.hx*nx + qk.hy*ny) / qk.h;

    // Lax-Wendroff local dissipation coefficient
    double vm = sqrt(grav*qm.h) + sqrt(qm.hx*qm.hx + qm.hy*qm.hy) / qm.h;
    double vk = sqrt(grav*qk.h) + sqrt(qk.hx*qk.hx + qk.hy*qk.hy) / qk.h;
    double a  = dt * std::max(vm*vm, vk*vk);

    // Helper values
    double scale = 0.5 * n_length;
    double gh2   = 0.5 * grav * (qm.h*qm.h + qk.h*qk.h);

    // Simple flux with dissipation for stability
    return QVar(scale * (wm*qm.h  + wk*qk.h)           - a * (qm.h  - qk.h),
                scale * (wm*qm.hx + wk*qk.hx + gh2*nx) - a * (qm.hx - qk.hx),
                scale * (wm*qm.hy + wk*qk.hy + gh2*ny) - a * (qm.hy - qk.hy));
  }
};

/** Node position function object for use in the SDLViewer. */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& n) {
    // positions of the nodes
    Point p = n.position();
    p.z = n.value().q.h;
    return p;
  }
};


template <typename MESH>
void output_debuging_info(MESH& m, double t) {
  for (auto it = m.triangle_begin(); it != m.triangle_end(); ++it){
    auto tri = *it;
    std::cerr << "Triangle " <<  tri.index() << " @" << t << "\n";
    std::cerr << "  Area " <<  tri.area() << "\n";
    std::cerr << "  Node positions (" <<  tri.node(0).position() << ")"
              << " (" <<  tri.node(1).position() << ")"
              << " (" <<  tri.node(2).position() << ")" << "\n";
    std::cerr << "  Triangle QVar [Q_bar, water column characteristics] h=" <<  tri.value().q_bar.h 
              << " hu=" << tri.value().q_bar.hx 
              << " hv=" << tri.value().q_bar.hy << "\n";

    // Edge Info
    for(int w = 0; w < 3; ++w){
      auto e = tri.edge(w);
      std::cerr << "  Edge " << w << " (" <<  tri.edge(w).node1().position() << ") (" 
                << tri.edge(w).node2().position() <<")\n";
      std::cerr << "    Normal (" <<  tri.normals_vector( e )  << ")\n";
      std::cerr << "    Adj Triangles" ;
        for (unsigned i=0; i<tri.edge(w).num_adj_triangles(); ++i)
          std::cerr << " " << tri.edge(w).triangle(i).index() ;
        std::cerr << "\n"; ;
      std::cerr << "    Opposite Tri QVar h=" << tri.edge(w).opposite_triangle(tri).value().q_bar.h;

      auto opp_tri = tri.edge(w).opposite_triangle(tri);
      // when the edge is a boundary edge, the triangle will be equial
      if(tri == opp_tri){
        // no opposite triangle... what to do?
      }
      std::cerr << "  hu=" << opp_tri.value().q_bar.hx 
                << " hv=" << opp_tri.value().q_bar.hy << "\n";
      

			if(e.num_adj_triangles() == 1)
				std::cerr << "    Boundary";
			else
				std::cerr << "    Edge";

			std::cerr << " flux h=" << tri.value().edge_fluxes[w].h
								<< " hu=" << tri.value().edge_fluxes[w].hx 
								<< " hv=" << tri.value().edge_fluxes[w].hy <<  "\n";
			
									
      /**
      // TODO: should we account for all fluxes in the fluxes vector?
      // I tried it, but the totals are different than Chris' values
      std::cerr << " flux h=" << e.value().total_flux().h
                << " hu=" << e.value().total_flux().hx 
                << " hv=" << e.value().total_flux().hy <<  "\n";**/
    }
    std::cerr << "\n";
  }
}

/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time.
 */
template <typename MESH, typename FLUX>
double hyperbolic_step(MESH& m, FLUX& f, double t, double dt) {

  // Step the finite volume model in time by dt.
  // Pseudocode:
  // Compute all fluxes. (before updating any triangle Q_bars)
  // For each triangle, update Q_bar using the fluxes as in Equation 8.
  // NOTE: Much like symp_euler_step, this may require TWO for-loops
  
  std::vector<QVar> sum_fluxes;// (m.num_triangles());
   
  //Iterate through all triangles and calculate new_q.
  for (auto it = m.triangle_begin(); it != m.triangle_end(); ++it){
    auto tri = *it;
    
    QVar sum_flux;
    tri.value().edge_fluxes.clear();
    assert(tri.value().edge_fluxes.size() == 0);

    //Iterate through each edge in this triangle
    for (auto i = 0; i <3; ++i) {
    	auto e = tri.edge(i);
    
    	//Boundary Conditions
      QVar adj_triangle_qbar;
      assert(e.num_adj_triangles() > 0 && e.num_adj_triangles() < 3);
      //std::cout << "num_adj_triangles " << e.num_adj_triangles() << std::endl;
      if (e.num_adj_triangles() == 1)
     		adj_triangle_qbar = QVar(tri.value().q_bar.h, 0, 0); 
      else
     		adj_triangle_qbar = e.opposite_triangle(tri).value().q_bar ; //adj_triangle_qbar = e.triangle(abs(i-1)).value().q_bar ; 	
     		
    	
    	QVar edge_flux = f( tri.normals_vector(e).x, 
                          tri.normals_vector(e).y, 
                          dt, 
                          tri.value().q_bar ,
                          adj_triangle_qbar);
      
      std::cout << "edge >>>>>> "<< edge_flux.h << " , " << edge_flux.hx << " , " << edge_flux.hy << std::endl;
      
      //for debugging purposes only
    	tri.value().edge_fluxes.push_back(edge_flux);	
    	sum_flux = sum_flux + edge_flux;

    }
    std::cout << "tri.value().edge_fluxes.size() >>>>>> "<< tri.value().edge_fluxes.size() << std::endl;
    
    tri.value().q_tmp = sum_flux;
    sum_fluxes.push_back(sum_flux);
    
    std::cout << "sum_flux: "<< sum_flux.h << " , " << sum_flux.hx << " , " << sum_flux.hy << std::endl;
    std::cout << std::endl;
    //std::cout << sum_flux.h << " , " << sum_flux.hx << " , " << sum_flux.hy << std::endl;
    //std::cout << tri.value().q_tmp.h << " , " << tri.value().q_tmp.hx << " , " << tri.value().q_tmp.hy << std::endl;
    //assert(sum_flux == tri.value().q_tmp);
  }
  
  output_debuging_info(m, t);
  //Iterate through all triangles and now update Qbar using fluxes calculated above.

  int l = 0;  
  for (auto it = m.triangle_begin(); it != m.triangle_end(); ++it){
    auto tri = *it;   
    // Equation 8
    if(l == 0){
      std::cout << "===========> "<< sum_fluxes[tri.index()].h << " , " << sum_fluxes[tri.index()].hx << " , " << sum_fluxes[tri.index()].hy << std::endl;
      auto prod = tri.value().q_tmp * (dt / tri.area());
      std::cout << "prod > "<< prod.h << " , " << prod.hx << " , " << prod.hy << std::endl;
      std::cout << "dt > " << dt << std::endl;
      std::cout << "dt/area > " << (dt / tri.area()) << std::endl;
      //std::cout << "           > "<< tri.value().q_tmp.h << " , " << tri.value().q_tmp.hx << " , " << tri.value().q_tmp.hy << std::endl;
    }

    ++l;
    tri.value().q_bar = tri.value().q_bar - (tri.value().q_tmp * (dt / tri.area()));
    //tri.value().q_bar = tri.value().q_bar - (sum_fluxes[tri.index()] * (dt / tri.area()));
    std::cout << "*** "<< tri.value().q_bar.h << " , " << tri.value().q_bar.hx << " , " << tri.value().q_bar.hy << std::endl;
      
  }
  
  return t + dt;
}
 
  // Translate the triangle-averaged values to node-averaged values
  // Implement Equation 9 from your pseudocode here


/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename MESH>
void post_process(MESH& m) {
	//Iterate through all the nodes in the mesh
	for (auto it = m.node_begin(); it != m.node_end(); ++it) {
    auto n = *it;  
    double area = 0;
    QVar q_k;
    
		//For each node, calculate the sum of q_k*tri_area for all adj triangles  
		for (auto t_it = n.triangle_begin(); t_it != n.triangle_end(); ++t_it) {
			auto t = *t_it;
			area = area + t.area();  
			q_k = q_k + (t.value().q_bar * t.area());
		}		
		n.value().q = q_k*(1.0 / area);
  }
}


/** Finds the edge with min length */
bool min_edge_length_func(edge_type e1, edge_type e2) { 
  return e1.length() < e2.length(); 
}

/** Define Initial Precondition */
struct InitialCondition {
  QVar operator()(const Point p){
    if (p.x < 0)
      return QVar(1.75, 0, 0); 
    
    return QVar(1, 0, 0);  
  }
};


int main(int argc, char* argv[])
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: shallow_water NODES_FILE TRIS_FILE\n";
    exit(1);
  }

  MeshType mesh;
  std::vector<node_type> mesh_node;

  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while (CS207::getline_parsed(tris_file, t)) {
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;


	// Set the initial values of the nodes and get the maximum height double
	double max_height = 0;
	auto init_cond = InitialCondition(); 
	for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it) { 
 	 	auto n = *it;
    n.value().q = init_cond(n.position());
  	max_height = std::max(max_height, n.value().q.h); 
	}
  std::cout << "- max height: " << max_height << std::endl;

  double min_edge_length = (*std::min_element(mesh.edge_begin(), mesh.edge_end(), min_edge_length_func)).length();
  std::cout << "- min length: " << min_edge_length << std::endl;

  // Set the initial values of the triangles to the average of their nodes
  for (auto it = mesh.triangle_begin(); it != mesh.triangle_end(); ++it) {
    auto t = *it;
    t.value().q_bar = (t.node(0).value().q +
                       t.node(1).value().q +
                       t.node(2).value().q) / 3.0;                                 
  }

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  /** Add nodes and edges to viewer */
  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                   CS207::DefaultColor(), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
  viewer.center_view();


  // CFL stability condition requires dt <= dx / max|velocity|
  // For the shallow water equations with u = v = 0 initial conditions
  //   we can compute the minimum edge length and maximum original water height
  //   to set the time-step
  // Compute the minimum edge length and maximum water height for computing dt
  double dt = 0.25 * min_edge_length / (sqrt(grav * max_height));
  std::cout << "- dt: " << dt << std::endl;


  /** Time variables */
  double t_start = 0;
  double t_end = 10;

  // Preconstruct a Flux functor
  EdgeFluxCalculator f;

  // Begin the time stepping
  for (double t = t_start; t < t_end; t += dt) {
    // Step forward in time with forward Euler
    hyperbolic_step(mesh, f, t, dt);

    // Update node values with triangle-averaged values
    post_process(mesh);

    // Update the viewer with new node positions
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                     CS207::DefaultColor(), NodePosition(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small meshes.
    // Feel free to remove them or tweak the constants.
    if (mesh.num_nodes() < 100)
      CS207::sleep(2);
  }

  return 0;
}
