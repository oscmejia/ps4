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
  
  /** Type value for nodes custom data */
  typedef N node_value_type;
  
  /** Type value for edge custom data */
  typedef E edge_value_type;
  
  
  /** Type of sizes. */
  typedef unsigned size_type;

  /** Type value for indexes */
  typedef unsigned idx_type;



  /** @struct InternalTriangle */
  struct InternalTriangle
  {
    idx_type node_idx1;
    idx_type node_idx2;
    idx_type node_idx3;
    idx_type center_idx;
    triangle_value_type value;
    

    InternalTriangle(idx_type node_idx1, idx_type node_idx2, idx_type node_idx3, triangle_value_type value) 
      : node_idx1(node_idx1), node_idx2(node_idx2), node_idx3(node_idx3), value(value) {
    }

  };
  /** Type of InternalTriangle */
  typedef InternalTriangle internal_triangle;

  
  /** Def of InternalNode */
  struct InternalNode
  {
    std::vector<idx_type> adj_triangles_;
  };
  /** Type of InternalNode */
  typedef InternalNode internal_node;


  /** Def of InternalEdge */
  struct InternalEdge
  {
    std::vector<idx_type> adj_triangles_;
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


  class NodeIncidentIterator;
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
        return std::tie(gn_.uid_, m_) < std::tie(x.gn_.uid_, x.m_);
      }

      /**
       * Returns node_incident_iterator pointing to the first element.
       * @return node_incident_iterator
       */
      node_incident_iterator triangle_begin() const {
        return NodeIncidentIterator(m_, gn_.value(), 0);
      }

      /**
       * Returns node_incident_iterator pointing to one element past the last valid element.
       * @return node_incident_iterator
       */
      node_incident_iterator triangle_end() const {
        void(); //return NodeIncidentIterator(m_, gn_.value(), m_->adj_n2t_[gn_.value()].size() );
      }
      
      // Return this node's value.
      node_value_type& value(){
      	return gn_.value();  //right? because gn_.value() only returns the value field within the internal_node 
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
  
  class Edge : public totally_ordered<Edge>  {
    public:

      /** Test whether this edge and @a x are equal.
       *
       * Equal edges are from the same mesh and have the same graph_edge.
       */
      bool operator==(const Edge& x) const {
        void(); //return std::tie(m_, ge_.node1_uid, ge_.node2_uid) == std::tie(x.m_, x.ge_.node1_uid, x.ge_.node2_uid);
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
        return std::tie(m_, ge_.node1_uid, ge_.node2_uid) < std::tie(x.m_, x.ge_.node1_uid, x.ge_.node2_uid);
      }

      triangle_type triangle(idx_type i) {
        assert(i < num_adj_triangles());
        // TODO: implement
      }

      size_type num_adj_triangles() {
        // TODO: implement
      }
      
      // Return this edge's value.
      edge_value_type& value(){
      	return ge_.value();  
      };

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
				void(); //return Edge(m_,m_->g_nodes_.has_edge(nodes_[(i+1)% 3], nodes_[(i+2)% 3]));
        
      }

      /** Accessing outward normal vectors of an edge of a triangle.*/
      Point normals_vector(const Edge& e) {
        // TODO: write implementation
      }

      double area() const {
        // http://www.mathopenref.com/coordtrianglearea.html
        Point a = nodes_[0].point;
        Point b = nodes_[1].point;
        Point c = nodes_[2].point;
        return std::abs( ( a.x*(b.y-c.y) + b.x*(c.y-a.y) + c.x*(a.y-b.y) ) / (2) );
      }

      /** Test whether this triangle and @a x are equal.
       *
       * Equal triangles are from the same mesh and have the same nodes and edges.
       */
      bool operator==(const Triangle& x) const {
        return std::tie(m_, nodes_[0], nodes_[1], nodes_[2]) == std::tie(x.m_, x.node(0), x.node(1), x.node(2));
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
        // TODO: can we just compare triangle idx ?
        return std::tie(m_, nodes_[0], nodes_[1], nodes_[2]) < std::tie(x.m_, x.node(0), x.node(1), x.node(2));
      }

      /**
       * Returns a reference to this triangle's value of type T.
       *
       * @return Object of type T by reference
       */
      triangle_value_type& value() {
        return m_->g_triangles_.node(idx_).value();
      }

      /**
       * Returns a reference to this trinagle's value of type T as a constant.
       *
       * @return Object of type T by reference as a constant.
       */
      const triangle_value_type& value() const {
        return m_->g_triangles_.node(idx_).value();
      }



      /**
       * Returns a Triangle adjacent to the edge idx provided.
       * @return triangle_type
       */
      triangle_type triangle(idx_type edge_idx) const {
        // TODO: implement
      }


      /**
       * Returns a Triangle adjacent to the edge provided.
       * @return triangle_type
       */
      triangle_type triangle(edge_type edge) const {
        // TODO: implement
      }


    private:
      
      Triangle(Mesh* m, Node node1, Node node2, Node node3, idx_type idx) 
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
    // TODO: clear mesh data structures
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

  /** Stores a graph node and returns a mesh node */
  Node add_node(const Point& position){
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

   
   	//Graph::Edge objects
    auto e1 = g_nodes_.add_edge(a.gn_, b.gn_);
    auto e2 = g_nodes_.add_edge(b.gn_, c.gn_);
    auto e3 = g_nodes_.add_edge(a.gn_, c.gn_);
    
    //Create Mesh::Edge objects for comparison with existing triangles to create adj_list below
    edge_type E1 = Edge(this,e1); 
    edge_type E2 = Edge(this,e2);
    edge_type E3 = Edge(this,e3);
    
    // TODO: we are storing an empty/invalid Point. other options? - I think that's fine, it's really a spaceholder for now.
    auto n = g_triangles_.add_node( Point(), InternalTriangle(a.gn_.index(), b.gn_.index(), c.gn_.index(), value));

    // TODO: Find adjacent triangles to this triangle, then call g_triangles_.add_edge() for each adjacent triangle found.
    
    //Iterate through all triangles. If there is a common edge, they are adjacent.  
		for (auto it = triangle_begin(); it != triangle_end(); ++it) {
			auto t = *it; 
			if (t.edge(0) == E1 || t.edge(0) == E2 || t.edge(0) == E3)
				g_triangles_.add_edge(n,g_triangles_.node(t.index()));
		}
	
    // TODO: Add idxs of triangles adjacent to the edges of this triangle 
    // TODO: Add idxs of triangles adjacent to the nodes of this triangle
    
    return Triangle(this, a, b, c, g_triangles_.num_nodes()-1);
  };









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
     * @Return Triangle object.
     */
    Triangle operator*() const {
      // TODO: construct triangle
      return Triangle();
    }

    /**
     * Incremental Operator for TriangleIterator.
     * Complexity: O(1).
     *
     * @Return TriangleIterator object by reference.
     */
    TriangleIterator& operator++() {
      //++idx_;
      return *this;
    }

    /**
     * Equality Operator for TriangleIterator.
     * Complexity: O(1).
     *
     * @Return bool, true if both TriangleIterator's are equial.
     */
    bool operator==(const TriangleIterator& other) const {
      void(); //return std::tie(m_, idx_) == std::tie(other.m_, other.idx_);
    }

   private:

    TriangleIterator(const GraphTriangleType* g)
        : gt_(const_cast<GraphTriangleType*>(g)){
    }

    /** Reference to the graph holding triangle nodes */
    GraphTriangleType* gt_;
    
    /** Node Iterator from g_triangles_*/
    //idx_type idx_;

    friend class Mesh;
  };

  /**
   * Return a triangle_iterator pointing to the begining
   * Complexity: O(1).
   *
   * @return TriangleIterator
   */
  triangle_iterator triangle_begin() const {
    return g_triangles_.node_begin(); // need to wrap around a triangle_iterator object...
  }

  /**
   * Return a triangle_iterator pointing to one pass the last valid position.
   * Complexity: O(1).
   *
   * @return TriangleIterator
   */
  triangle_iterator triangle_end() const {
    void(); //return triangle_iterator( this, g_triangles_.size() );
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
      // TODO: construct node
      return Node();
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
      // TODO: construct edge
      return Edge();
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






 private:

  /** Graph that holds nodes and edges*/
  GraphType g_nodes_;

  /** Graph that holds nodes that represent the center of each triangle and holds internal_triangle objects*/
  GraphTriangleType g_triangles_;


};

#endif