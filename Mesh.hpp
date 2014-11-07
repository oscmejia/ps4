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
template <typename T>
class Mesh {

 public:

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** Type of this mesh. */
  typedef Mesh mesh_type;

  /** Type value for triangles custom data */
  typedef T triangle_value_type;

  /** Type of indexes and sizes. */
  typedef unsigned size_type;

  /** Type value for idexes */
  typedef unsigned idx_type;

  /** Define types from Graph */
  typedef Graph<int,int> GraphType;
  typedef typename GraphType::Node GNode;
  typedef typename GraphType::Edge GEdge;
  
   
  /** @struct InternalTriangle */
  struc InternalTriangle
  {
    idx_type node_idx1;
    idx_type node_idx2;
    idx_type node_idx3;
    idx_type center_idx;
    triangle_value_type value;
    

    InternalTriangle(idx_type node_idx1, idx_type node_idx2, idx_type node_idx3, size_type center_idx, triangle_value_type value) 
      : node_idx1(node_idx1), node_idx2(node_idx2), node_idx3(node_idx3), center_idx(center_idx), value(value) {
    }

  };

  /** Type of InternalTriangle */
  typedef InternalTriangle internal_triangle;

  /** Type of triangle iterators, which iterate over all mesh triangles. */
  class TriangleIterator;
  /** Synonym for TriangleIterator */
  typedef TriangleIterator triangle_iterator;


  /** Mesh Node*/
  class Node;
  /** Type of Mesh Node*/
  typedef Node node_type;
  

  /** Mesh Edge*/
  class Edge;
  /** Type of Mesh Edge*/
  typedef Edge edge_type;


  /** Type of incident iterators, which iterate incident triangles to an edge. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Construct an empty Mesh. */
  Mesh(){
  }

  /** Default destructor */
  ~Mesh() = default;


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

  /** Accessing outward normal vectors of an edge of a triangle.*/
  // TODO: validate signature
  std:vector<double> normals_vector(const Triangle& t, const& Edge e) {

  }


  ////////////////////
  // MESH NODES     //
  ////////////////////
  
  class Node : public GNode<Node>  {

  };




  ////////////////////
  // MESH EDGES     //
  ////////////////////
  
  class Edges : public GEdge<Edges>  {
    public:

      /**
       * Returns incident_iterator poiting to the first element.
       * @return incident_iterator
       */
      incident_iterator triangle_begin() const {
        return IncidentIterator(g_, uid_, 0);
      }

      /**
       * Returns incident_iterator poiting to one elemnt past the last valid element.
       * @return incident_iterator
       */
      incident_iterator triangle_end() const {
        return IncidentIterator(g_, uid_, g_->i2u_edges_[index()].size());
      }


      void add adj_triangle(idx_type idx) {

      }

    private:
      Mesh m_;

      Edge(const Mesh* m, size_type uid) : m_(const_cast<Mesh*>(m)), uid_(uid) {
        assert(m_ != nullptr);
      }


      std::vector<idx_type> adj_triangles_
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
      /** Construct an invalid Edge. */
      Triangle() {
      }

      /** Return node 1 of this Triangle */
      Node node1() const {
        return node1_;
      }

      /** Return node 2 of this Triangle */
      Node node2() const {
        return node2_;
      }

      /** Return node 3 of this Triangle */
      Node node3() const {
        return node3_;
      }

      Edge edge1() const {
        return 
      }

      Edge edge2() const {
        
      }

      Edge edge3() const {
        
      }


      double area() const {
        // http://www.mathopenref.com/coordtrianglearea.html
        Point a = node1_.point;
        Point b = node2_.point;
        Point c = node3_.point;
        return std::abs( ( a.x*(b.y-c.y) + b.x*(c.y-a.y) + c.x*(a.y-b.y) ) / (2) );
      }


      /** Test whether this triangle and @a x are equal.
       *
       * Equal triangles are from the same mesh and have the same nodes and edges.
       */
      bool operator==(const Triangle& x) const {
        return std::tie(m_, uid1_, uid2_, uid3_) == std::tie(x.m_, x.uid1_, x.uid2_, x.uid3_);
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
        // TODO: can we just compare idx ?
        return std::tie(m_, node1_, node2_, node3_) < std::tie(x.m_, x.node1_, x.node2_, x.node3_);
      }

      double length () const {
        // TODO
        //return norm(node1().position() - node2().position());
      }

      /**
       * Returns a reference to this triangle's value of type T.
       *
       * @return Object of type T by reference
       */
      triangle_value_type& value() {
        return  m_->internal_triangles_[idx_].value;
      }

      /**
       * Returns a reference to this trinagle's value of type T as a constant.
       *
       * @return Object of type T by reference as a constant.
       */
      const triangle_value_type& value() const {
        return  m_->internal_triangles_[idx_].value;
      }


    private:
      
      Triangle(Mesh* m, Node node1, Node node2, Node node3, idx_type idx) 
        : m_(const_cast<Mesh*>(m)), node1_(node1), node2_(node1), node3_(node1), idx_(idx) {
      }

      Mesh* m_;
      Node node1_;
      Node node2_;
      Node node3_;
      idx_type idx_;

      friend class Mesh;

  };

  /** Add an triangle to the graph, or return the current triangle if it already exists.
   * @pre @a a, @a b and @a c are distinct points that confirm a triangle
   * @return an Triangle object e with t.node1() == @a a, t.node2() == @a b and t.node3() == @a c
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Triangle add_triangle(const Point& a, const Point& b, const Point& c, const triangle_value_type& value = triangle_value_type()) {
    assert(a != b && b != c && a != c);
    auto n1 = g_nodes_.add_node(a);
    auto n2 = g_nodes_.add_node(b);
    auto n3 = g_nodes_.add_node(c);

    // TODO: do we really need to add edges?
    g_nodes_.add_edge(a, b);
    g_nodes_.add_edge(b, c);
    g_nodes_.add_edge(a, c);


    // TODO: add triangles to g_triangles. We need to find adjacent triangles first
    // TODO: add center node to g_triangles_ 
    // TODO: replace 0 in internal_triangles and Triangle constructor
    
    internal_triangles_.push_back( InternalTriangle(n1.index(), n2.index(), n3.index(), 0, value) );

    return Triangle(this, n1, n2, n3, 0);
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
   * Return a node_iterator pointing to the begining
   * Complexity: O(1).
   *
   * @return TriangleIterator
   */
  triangle_iterator triangle_begin() const {
    return triangle_iterator( this, 0 );
  }

  /**
   * Return a node_iterator pointing to one pass the last valid position.
   * Complexity: O(1).
   *
   * @return TriangleIterator
   */
  triangle_iterator triangle_end() const {
    return triangle_iterator( this, g_triangles_.size() );
  }






 private:

  /** Graph that holds nodes*/
  GraphType g_nodes_;

  /** Graph that holds nodes that represent the center of each triangle*/
  GraphType g_triangles_;

  /** Stores all triangle objects, indexed by triangle idxs */
  std::vector<internal_triangle> internal_triangles_;

  
  

};

#endif