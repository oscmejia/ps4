#ifndef CS207_MESH_HPP
#define CS207_MESH_HPP

#include "Graph.hpp"
/**
 * @brief A Mesh
 */


/** @class Mesh
 * @brief A template for a Mesh.
 *
 * RI(mesh):
 *
 * Users can add and retrieve nodes, edges and triangles.
 */
template <typename T, typename N, typename E>
class Mesh {

 public:

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** Type of this mesh. */
  typedef Mesh mesh_type;

  /** Type value for triangles custom data */
  typedef T triangle_value_type;
  
  /** Type value for nodes user data */
  typedef N node_value_type;
  
  /** Type value for edge user data */
  typedef E edge_value_type;
  
  /** Type of sizes. */
  typedef unsigned size_type;

  /** Type value for indexes */
  typedef unsigned idx_type;


  /** @struct InternalTriangle */
  struct InternalTriangle {
    idx_type node_idx1;
    idx_type node_idx2;
    idx_type node_idx3;
    triangle_value_type user_value;
    
    InternalTriangle(idx_type node_idx1, idx_type node_idx2, idx_type node_idx3, triangle_value_type value) 
      : node_idx1(node_idx1), node_idx2(node_idx2), node_idx3(node_idx3), user_value(value) {
    }
  };
  /** Type of InternalTriangle */
  typedef InternalTriangle internal_triangle;

  
  /** Def of InternalNode */
  struct InternalNode {
    std::vector<idx_type> adj_triangles_;
    node_value_type user_value;
    
    /** Constructor with user value */
    InternalNode(node_value_type value) 
      : user_value(value) {
    }
    
    //Default Constructor
    InternalNode() {
    }
  };
  /** Type of InternalNode */
  typedef InternalNode internal_node;


  /** Def of InternalEdge */
  struct InternalEdge {
    std::vector<idx_type> adj_triangles_;
    edge_value_type user_value;
    
  	/** Constructor with user value */
    InternalEdge(node_value_type value) 
      : user_value(value) {
    }
    
    /** Default Constructor */
    InternalEdge() {
    }
  };
  /** Type of InternalEdge */
  typedef InternalEdge internal_edge;


  /** Define types from Graph container for nodes/edges */
  typedef Graph<internal_node, internal_edge> GraphType;
  typedef typename GraphType::Node graph_node;
  typedef typename GraphType::Edge graph_edge;
  typedef typename GraphType::NodeIterator graph_node_iterator;
  typedef typename GraphType::EdgeIterator graph_edge_iterator;
  /** Define types from Graph container for triangles */
  typedef Graph<internal_triangle, int> GraphTriangleType;


  /** Type of node incident iterators, which iterate over all triangles adjacent to a node. */
  class NodeIncidentIterator;
  /** Synonym for NodeIncidentIterator */
  typedef NodeIncidentIterator node_incident_iterator;


  /** Type of triangle iterators, which iterate over all mesh triangles. */
  class TriangleIterator;
  /** Synonym for TriangleIterator */
  typedef TriangleIterator triangle_iterator;


  /** Type of node iterators, which iterate over all mesh nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;


  /** Type of edge iterators, which iterate over all mesh edges. */
  class EdgeIterator;
  /** Synonym for NodeIterator */
  typedef EdgeIterator edge_iterator;


  /** Mesh Triangle*/
  class Triangle;
  /** Type of Mesh Triangle*/
  typedef Triangle triangle_type;


  /** Mesh Node*/
  class Node;
  /** Type of Mesh Node*/
  typedef Node node_type;
  

  /** Mesh Edge*/
  class Edge;
  /** Type of Mesh Edge*/
  typedef Edge edge_type;


  /** Type of incident iterators, which iterate incident triangles to an edge. */
  class EdgeIncidentIterator;
  /** Synonym for IncidentIterator */
  typedef EdgeIncidentIterator edge_incident_iterator;


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Construct an empty Mesh. */
  Mesh(){
  }

  /** Default destructor */
  ~Mesh() = default;


  ////////////////////
  // MESH NODES     //
  ////////////////////
  
  /** @class Mesh::Node
   * @brief Class representing mesh's nodes.
   *
   */
  class Node : public totally_ordered<Node>  {
    public:

      /** Test whether this nodes and @a x are equal.
       *
       * Equal nodes are from the same mesh and have the same graph_edge.
       */
      bool operator==(const Node& x) const {
        return std::tie(gn_, m_) == std::tie(x.gn_, x.m_);
      }

      /** Test whether this node is less than @a x in the global order.
       *
       * This ordering function is useful for STL containers such as
       * std::map<>. It need not have any geometric meaning.
       *
       * The node ordering relation must obey trichotomy: For any two nodes x
       * and y, exactly one of x == y, x < y, and y < x is true.
       */
      bool operator<(const Node& x) const {
        return std::tie(gn_, m_) < std::tie(x.gn_, x.m_);
      }
      
      //Returns position of current Mesh::Node
      Point position() const {
      	return gn_.position();
      }

      /**
       * Returns node_incident_iterator pointing to the first element.
       * @return node_incident_iterator
       */
      node_incident_iterator triangle_begin() const {
        return NodeIncidentIterator(m_, gn_.index(), 0);
      }

      /**
       * Returns node_incident_iterator pointing to one element past the last valid element.
       * @return node_incident_iterator
       */
      node_incident_iterator triangle_end() const {
        return NodeIncidentIterator(m_, gn_.index(), gn_.value().adj_triangles_.size() );
      }
      
      
      // Return this node's value.
      node_value_type& value() {
        return gn_.value().user_value;  
      };

      // Return this node's value.
      node_value_type& value() const {
        auto value =  gn_.value();
        auto uv = value.user_value;
        return uv;    
      };
      


    private:
      Mesh* m_;
      graph_node gn_;

      Node(const Mesh* m, graph_node gn) : m_(const_cast<Mesh*>(m)), gn_(gn) {
        assert(m_ != nullptr);
      }
      
      friend class Mesh;
  };


  ////////////////////
  // MESH EDGES     //
  ////////////////////
  
  /** @class Mesh::Edge
   * @brief Class representing mesh's edges.
   *
   */
  class Edge : public totally_ordered<Edge>  {
    public:

      /** Test whether this edge and @a x are equal.
       *
       * Equal edges are from the same mesh and have the same graph_edge.
       */
      bool operator==(const Edge& x) const {
        return std::tie(m_, ge_) == std::tie(x.m_, x.ge_);
      }

      /** Test whether this edge is less than @a x in the global order.
       *
       * This ordering function is useful for STL containers such as
       * std::map<>. It need not have any geometric meaning.
       *
       * The edge ordering relation must obey trichotomy: For any two edges x
       * and y, exactly one of x == y, x < y, and y < x is true.
       */
      bool operator<(const Edge& x) const {
        return std::tie(m_, ge_) < std::tie(x.m_, x.ge_);
      }

      triangle_type triangle(idx_type i) {
        assert(i < num_adj_triangles());   
        return m_->triangle(ge_.value().adj_triangles_[i]);
      }

      size_type num_adj_triangles() {
        return ge_.value().adj_triangles_.size();
      }
      
      // Return this edge's value.
      edge_value_type& value(){
      	return ge_.value().user_value;  
      };
      
      node_type node1() {
      	return Node(m_,ge_.node1());
    	};

		  node_type node2() {
      	return Node(m_,ge_.node2());
      };

      double length() const {
        return ge_.length();
      }

    private:

      Mesh* m_;

      graph_edge ge_;

      Edge(const Mesh* m, graph_edge ge) 
        : m_(const_cast<Mesh*>(m)), ge_(ge) {

        assert(m_ != nullptr);
      }
      
      friend class Mesh;
  };


  ////////////////////
  // MESH TRIANGLES //
  ////////////////////

  /** @class Mesh::Triangle
   * @brief Class representing mesh's triangles.
   *
   */
  class Triangle : private totally_ordered<Triangle>  {
    public:

      /** Construct an invalid Triangle. */
      Triangle() {
      }
      
      /** Return this triangles's index, a number in the range [0, num_triangles). */
      size_type index() const {
        return idx_;
      }

      /** Return node at position i in this Triangle */
      node_type node(idx_type i) const {
        assert(i >= 0 && i < 3);
        return nodes_[i];
      }

			//Return edge i, defined as the one opposite node i, i.e composed of the other two nodes in the triangle
      edge_type edge(idx_type i) const {
        assert(i >= 0 && i < 3);

        //TODO: returning the nodes in this way so is easy to compare/debug
				if(i == 0)
          return Edge(m_, m_->g_nodes_.add_edge(nodes_[0].gn_, nodes_[1].gn_) );
        
        if(i == 1)
          return Edge(m_, m_->g_nodes_.add_edge(nodes_[1].gn_, nodes_[2].gn_) );
        
        if(i == 2)
          return Edge(m_, m_->g_nodes_.add_edge(nodes_[2].gn_, nodes_[0].gn_) );
        
        //return Edge(m_, m_->g_nodes_.add_edge(nodes_[(i+1)% 3].gn_,nodes_[(i+2)% 3].gn_)); 
      }

      /** Accessing outward normal vectors of an edge of a triangle.*/
      Point normals_vector(Edge& e) {
      	Point p1 = e.node1().position();
      	Point p2 = e.node2().position();
      	
      	Point edge = p1-p2;
      	
        //TODO: This prevents the z negative zeros, but nothing else. revert after all is working.
        return Point( edge.y - edge.z,
                      -edge.x,
                      0);

      	//return cross(edge, Point(0,0,1));
    	}

      /** Calculates the area of a triangle */
      double area() const {
        Point a = nodes_[0].position();
        Point b = nodes_[1].position();
        Point c = nodes_[2].position();
        return std::abs( ( a.x*(b.y-c.y) + b.x*(c.y-a.y) + c.x*(a.y-b.y) ) / (2) );
      }

      /** Test whether this triangle and @a x are equal.
       *
       * Equal triangles are from the same mesh and have the same nodes and edges.
       */
      bool operator==(const Triangle& x) const {
        //return std::tie(m_, nodes_[0], nodes_[1], nodes_[2]) == std::tie(x.m_, x.node(0), x.node(1), x.node(2));
        return std::tie(m_,idx_) == std::tie(x.m_,x.idx_);
      }

      /** Test whether this triangle is less than @a x in the global order.
       *
       * This ordering function is useful for STL containers such as
       * std::map<>. It need not have any geometric meaning.
       *
       * The edge ordering relation must obey trichotomy: For any two edges x
       * and y, exactly one of x == y, x < y, and y < x is true.
       */
      bool operator<(const Triangle& x) const {
        return std::tie(idx_,m_) < std::tie(x.idx_,m_);
      }

      /**
       * Returns a reference to this triangle's user value of type T.
       *
       * @return Object of type T by reference
       */
      triangle_value_type& value() {
        return m_->g_triangles_.node(idx_).value().user_value;
      }

      /**
       * Returns a reference to this triangle's value of type T as a constant.
       *
       * @return Object of type T by reference as a constant.
       */
      const triangle_value_type& value() const {
        return m_->g_triangles_.node(idx_).value().user_value;
      }

    private:
      
      Triangle(const Mesh* m, Node node1, Node node2, Node node3, idx_type idx) 
        : m_(const_cast<Mesh*>(m)), idx_(idx) {
        
        nodes_.push_back(node1);
        nodes_.push_back(node2);
        nodes_.push_back(node3);
      }

      Mesh* m_;
      std::vector<Node> nodes_;
      idx_type idx_;

      friend class Mesh;
  };



  /////////////
  // GENERAL //
  /////////////

  /** Remove all nodes, edges and triangles from this graph.
   * @post num_nodes() == 0 && num_edges() == 0 && num_triangles() == 0
   *
   * Invalidates all outstanding Node, Edge and Triangle objects.
   */
  void clear() {
    g_nodes_.clear();
    g_triangles_.clear();
  }


  /** Returns the number of triangles */
  size_type num_triangles() const {
    return g_triangles_.num_nodes();
  }

  /** Returns the number of nodes */
  size_type num_nodes() const {
    return g_nodes_.num_nodes();
  }

  /** Returns the number of edges */
  size_type num_edges() const {
    return g_nodes_.num_edges();
  }
  
  
	/** Return triangle at position i in this Mesh */
	triangle_type triangle(idx_type i) const {
		assert(i >= 0 && i < num_triangles());
	
		auto n1_idx= g_triangles_.node(i).value().node_idx1;
		auto n2_idx = g_triangles_.node(i).value().node_idx2;
		auto n3_idx = g_triangles_.node(i).value().node_idx3;
	
		return Triangle(this,
      Node(this,g_nodes_.node(n1_idx)),
      Node(this,g_nodes_.node(n2_idx)),
      Node(this,g_nodes_.node(n3_idx)),
      i);
	}


  /** Stores a graph node and returns a mesh node */
  Node add_node(const Point& position) {
    return Node(this, g_nodes_.add_node(position));
  }

  /** Add an triangle to the graph, or return the current triangle if it already exists.
   * @pre @a a, @a b and @a c are distinct nodes that form a triangle
   * @return an Triangle object e with t.node1() == @a a, t.node2() == @a b and t.node3() == @a c
   * @post has_edge(@a a, @a b) == true, has_edge(@a,@c) == true, has_edge(@b,@c) == true;
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Triangle add_triangle(const Node& a, const Node& b, const Node& c, 
    const triangle_value_type& value = triangle_value_type()) {
    
    assert(a != b && b != c && a != c);
    assert( g_nodes_.has_node(a.gn_) );
    assert( g_nodes_.has_node(b.gn_) );
    assert( g_nodes_.has_node(c.gn_) );

    
    //TODO: If triangle has already been added, return the Triangle
     

    //Create Mesh::Edge objects for comparison with existing triangles to create adj_list below
    edge_type e1 = Edge(this,g_nodes_.add_edge(a.gn_, b.gn_)); 
    edge_type e2 = Edge(this,g_nodes_.add_edge(b.gn_, c.gn_));
    edge_type e3 = Edge(this,g_nodes_.add_edge(c.gn_, a.gn_));
    
    //Add node to triangle graph
    auto n = g_triangles_.add_node( Point(), InternalTriangle(a.gn_.index(), b.gn_.index(), c.gn_.index(), value));

	 //Iterate through all triangles. If they are adjacent (have at least one common edge), add edge between triangle nodes.  
		for (auto it = triangle_begin(); it != triangle_end(); ++it) {
			auto t = *it; 
			//Cycle through t's 3 edges and compare with the current triangle's edges. 
			for (int j = 0; j<3; ++j){
				if (t.edge(j) == e1 || t.edge(j) == e2 || t.edge(j) == e3){ 
 					if (n != g_triangles_.node(t.index()))
						g_triangles_.add_edge(n,g_triangles_.node(t.index()));
					}
			}
		}
		
		idx_type tri_idx = g_triangles_.num_nodes()-1;
		
		//Adding new triangle to the adj_list of the three new edges (make sure this triangle has not been added before!)
		e1.ge_.value().adj_triangles_.push_back(tri_idx);
		e2.ge_.value().adj_triangles_.push_back(tri_idx);
		e3.ge_.value().adj_triangles_.push_back(tri_idx);
	
		//TODO : make code below work!
		
		//Adding new triangle idx to the adj_list of the three nodes (make sure this triangle has not been added before!)	
		  
		 // std::cerr<< "Size is " << a.gn_.value().adj_triangles_.push_back(1) ; //size();	
		auto v = a.gn_.value().adj_triangles_;
		v.push_back(tri_idx);
		
	  auto v2= b.gn_.value().adj_triangles_;
	  v2.push_back(tri_idx);
	  
		auto v3=c.gn_.value().adj_triangles_;
		v3.push_back(tri_idx);
			 
    return Triangle(this, a, b, c, tri_idx);
  }

  ///////////////
  // Iterators //
  ///////////////

  /** @class Mesh::TriangleIterator
   * @brief Iterator class for triangles. A forward iterator. */
  class TriangleIterator : private totally_ordered<TriangleIterator> {

   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

     /** Construct an invalid TriangleIterator. */
    TriangleIterator() {
    }

    /**
     * Reference operator for TriangleIterator.
     * Complexity: O(1).
     *
     * @Return Node object.
     */
    Triangle operator*() const {
      auto n1_idx= m_->g_triangles_.node(idx_).value().node_idx1;
      auto n2_idx = m_->g_triangles_.node(idx_).value().node_idx2;
      auto n3_idx = m_->g_triangles_.node(idx_).value().node_idx3;
      
      return Triangle(m_,
                      Node(m_,m_->g_nodes_.node(n1_idx)),
                      Node(m_,m_->g_nodes_.node(n2_idx)),
                      Node(m_,m_->g_nodes_.node(n3_idx)),
                      idx_);
    }

    /**
     * Incremental Operator for TriangleIterator.
     * Complexity: O(1).
     *
     * @Return TriangleIterator object by reference.
     */
    TriangleIterator& operator++() {
      ++idx_;
      return *this;
    }

    /**
     * Equialy Operator for TriangleIterator.
     * Complexity: O(1).
     *
     * @Return bool, true if both TriangleIterator's are equial.
     */
    bool operator==(const TriangleIterator& other) const {
      return std::tie(m_, idx_) == std::tie(other.m_, other.idx_);
    }

   private:

    TriangleIterator(const Mesh* m, idx_type idx)
      : m_(const_cast<Mesh*>(m)), idx_(idx) {
    }

    /** Reference to the mesh */
    Mesh* m_;
    /** Triangle idx */
    idx_type idx_;

    friend class Mesh;
  };

 /**
   * Return a triangle_iterator pointing to the begining
   * Complexity: O(1).
   *
   * @return TriangleIterator
   */
  triangle_iterator triangle_begin() const {
    return triangle_iterator( this, 0 );
  }

  /**
   * Return a triangle_iterator pointing to one pass the last valid position.
   * Complexity: O(1).
   *
   * @return TriangleIterator
   */
  triangle_iterator triangle_end() const {
    return triangle_iterator( this, g_triangles_.size() );
  }





  /** @class Mesh::NodeIterator
   * @brief Iterator class for Mesh::Node. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {

   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /**
     * Reference operator for NodeIterator.
     * Complexity: O(1).
     *
     * @Return Node object.
     */
    Node operator*() const {
      return Node(m_,*it_);
    }

    /**
     * Incremental Operator for NodeIterator.
     * Complexity: O(1).
     *
     * @Return NodeIterator object by reference.
     */
    NodeIterator& operator++() {
      ++it_;
      return *this;
    }

    /**
     * Equialy Operator for NodeIterator.
     * Complexity: O(1).
     *
     * @Return bool, true if both NodeIterators are equial.
     */
    bool operator==(const NodeIterator& other) const {
      return std::tie(m_, it_) == std::tie(other.m_, other.it_);
    }

   private:

   NodeIterator(const Mesh* m, graph_node_iterator it)
        : m_(const_cast<Mesh*>(m)), it_(it) {
    }

    /** Reference to the mesh */
    Mesh* m_;
    /** Graph::NodeIterator */
    graph_node_iterator it_;

    friend class Mesh;
  };

  /**
   * Return a node_iterator pointing to the begining
   * Complexity: O(1).
   *
   * @return NodeIterator
   */
  node_iterator node_begin() const {
    return node_iterator( this, g_nodes_.node_begin() );
  }

  /**
   * Return a node_iterator pointing to one pass the last valid position.
   * Complexity: O(1).
   *
   * @return NodeIterator
   */
  node_iterator node_end() const {
    return node_iterator( this, g_nodes_.node_end() );
  }






  /** @class Mesh::EdgeIterator
   * @brief Iterator class for Mesh::Edge. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {

   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /**
     * Reference operator for EdgeIterator.
     * Complexity: O(1).
     *
     * @Return Edge object.
     */
    Edge operator*() const {
      return Edge(m_,*it_);
    }

    /**
     * Incremental Operator for EdgeIterator.
     * Complexity: O(1).
     *
     * @Return EdgeIterator object by reference.
     */
    EdgeIterator& operator++() {
      ++it_;
      return *this;
    }

    /**
     * Equialy Operator for EdgeIterator.
     * Complexity: O(1).
     *
     * @Return bool, true if both EdgeIterators are equial.
     */
    bool operator==(const EdgeIterator& other) const {
      return std::tie(m_, it_) == std::tie(other.m_, other.it_);
    }

   private:

   EdgeIterator(const Mesh* m, graph_edge_iterator it)
        : m_(const_cast<Mesh*>(m)), it_(it) {
    }

    /** Reference to the mesh */
    Mesh* m_;
    /** Graph::EdgeIterator */
    graph_edge_iterator it_;

    friend class Mesh;
  };

  /**
   * Return a edge_iterator pointing to the begining
   * Complexity: O(1).
   *
   * @return NodeIterator
   */
  edge_iterator edge_begin() const {
    return edge_iterator( this, g_nodes_.edge_begin() );
  }

  /**
   * Return a edge_iterator pointing to one pass the last valid position.
   * Complexity: O(1).
   *
   * @return TriangleIterator
   */
  edge_iterator edge_end() const {
    return edge_iterator( this, g_nodes_.edge_end() );
  }


  /** @class Mesh::NodeIncidentIterator
   * @brief Iterator class for adjacnet triangles to a node. A forward iterator. */
  class NodeIncidentIterator : private totally_ordered<NodeIncidentIterator> {

   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

     /** Construct an invalid NodeIncidentIterator. */
    NodeIncidentIterator() {
    }

    /**
     * Reference operator for NodeIncidentIterator.
     * Complexity: O(1).
     *
     * @Return Triangle object.
     */
    Triangle operator*() const {

      // InternalNode
      auto this_node = m_->g_nodes_.node(node_idx_);

      // InternalTriangle index
      idx_type this_tri_idx = this_node.value().adj_triangles_[adj_tri_idx_];

      // InternalTriangle
      auto this_tri = m_->g_triangles_.node(this_tri_idx);

      return Triangle(m_,
                      Node(m_, m_->g_nodes_.node(this_tri.value().node_idx1)),
                      Node(m_, m_->g_nodes_.node(this_tri.value().node_idx2)),
                      Node(m_, m_->g_nodes_.node(this_tri.value().node_idx3)),
                      this_tri_idx);
    }

    /**
     * Incremental Operator for NodeIncidentIterator.
     * Complexity: O(1).
     *
     * @Return NodeIncidentIterator object by reference.
     */
    NodeIncidentIterator& operator++() {
      ++adj_tri_idx_;
      return *this;
    }

    /**
     * Equialy Operator for NodeIncidentIterator.
     * Complexity: O(1).
     *
     * @Return bool, true if both NodeIncidentIterator's are equial.
     */
    bool operator==(const NodeIncidentIterator& other) const {
      return std::tie(m_, node_idx_, adj_tri_idx_) == std::tie(other.m_, other.node_idx_, other.adj_tri_idx_);
    }

   private:

    NodeIncidentIterator(const Mesh* m, idx_type node_idx, idx_type adj_tri_idx)
        : m_(const_cast<Mesh*>(m)), node_idx_(node_idx), adj_tri_idx_(adj_tri_idx) {
    }

    /** Reference to the mesh */
    Mesh* m_;
    /** Node idx */
    idx_type node_idx_;
    /** Triangle idx */
    idx_type adj_tri_idx_;

    friend class Mesh;
  };


 private:

  /** Graph that holds nodes and edges*/
  GraphType g_nodes_;

  /** Graph that holds nodes that represent the center of each triangle and holds internal_triangle objects*/
  GraphTriangleType g_triangles_;


};

#endif