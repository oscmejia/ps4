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
template <typename V, typename E, typename T>
class Mesh : Graph<V, E> {

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
  typedef Graph<V,E> GraphType;
  typedef typename Graph<V,E>::Node Node;
  typedef typename Graph<V,E>::Edge Edge;
  typedef typename Graph<V,E>::node_type node_type;
  typedef typename Graph<V,E>::edge_type edge_type;
  
  /** Type of triangle iterators, which iterate over all mesh triangles. */
  class TriangleIterator;
  /** Synonym for TriangleIterator */
  typedef TriangleIterator triangle_iterator;

  /** Type of node2triangle iterators, which iterate over all triangles adjacent to a node. */
  class Node2TriangleIterator;
  /** Synonym for Node2TriangleIterator */
  typedef Node2TriangleIterator node2triangle_iterator;

  /** Type of edge2triangle iterators, which iterate over all triangles adjacent to an edge. */
  class Edge2TriangleIterator;
  /** Synonym for Edge2TriangleIterator */
  typedef Edge2TriangleIterator edge2triangle_iterator;


  /** @struct InternalTriangle */
  struct InternalTriangle
  {
    size_type uid1;
    size_type uid2;
    size_type uid3;
    idx_type idx;
    triangle_value_type value;

    InternalTriangle(size_type uid1, size_type uid2, size_type uid3, size_type idx, triangle_value_type value) 
      : uid2(uid1), uid2(uid2), uid3(uid3), idx(idx), value(value) {
    }
  };

  /** Type of InternalTriangle */
  typedef InternalTriangle internal_triangle;

  /** Inner vector of idxs for nodes to triangles */
  typedef std::vector<idx_type> idx_list_type;

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
    std::cout << "Mesh.clean()" << std::endl;
    // call g.clear()
  }


  /** Returns the number of triangles */
  size_type num_triangles() const {
    return num_triangles_();
  }


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
        return Node(g_, node1_uid);
      }

      /** Return node 2 of this Triangle */
      Node node2() const {
        return Node(g_, node2_uid);
      }

      /** Return node 3 of this Triangle */
      Node node3() const {
        return Node(g_, node3_uid);
      }

      /**
       * Return an edge for 2 triangle nodes 
       * @param  a node
       * @param  b node
       * @return   Edge
       */
      Edge edge(Node a, Node b) const {
        assert(i2u_nodes_[a.index()] == node1_uid || 
          i2u_nodes_[a.index()] == node2_uid || 
          i2u_nodes_[a.index()] == node3_uid);
        assert(i2u_nodes_[b.index()] == node1_uid || 
          i2u_nodes_[b.index()] == node2_uid || 
          i2u_nodes_[b.index()] == node3_uid);

        return Edge(this, i2u_nodes_[a.index()], i2u_nodes_[b.index()]);
      }


      double area() const {

      }


      /** Test whether this triangle and @a x are equal.
       *
       * Equal triangles are from the same mesh and have the same nodes and edges.
       */
      bool operator==(const Triangle& x) const {
        return std::tie(m_, node1_uid, node2_uid, node3_uid) == std::tie(x.m_, x.node1_uid, x.node2_uid, x.node3_uid);
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
        return std::tie(m_, node1_uid, node2_uid, node3_uid) < std::tie(x.m_, x.node1_uid, x.node2_uid, x.node3_uid);
      }

      double length () const {
        //return norm(node1().position() - node2().position());
      }

      /**
       * Returns a reference to this triangle's value of type T.
       *
       * @return Object of type T by reference
       */
      triangle_value_type& value() {
        return  internal_triangles_[t_uid].value;
      }

      /**
       * Returns a reference to this trinagle's value of type T as a constant.
       *
       * @return Object of type T by reference as a constant.
       */
      const triangle_value_type& value() const {
        return  internal_triangles_[t_uid].value;
      }


    private:
      
      Triangle(Mesh* m, size_type t_uid, size_type uid1, size_type uid2, size_type uid3) 
        : m_(const_cast<Mesh*>(m)), t_uid_(t_uid) uid1_(uid1), uid2_(uid2), uid3_(uid3) {
      }

      Mesh* m_;
      size_type uid1_;
      size_type uid2_;
      size_type uid3_;
      size_type t_uid_;

      friend class Mesh;

  };

  /** Add an triangle to the graph, or return the current triangle if it already exists.
   * @pre @a a, @a b and @a c are distinct valid nodes of the graph
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
  Triangle add_triangle(const Node& a, const Node& b, const Node& c, const triangle_value_type& value = triangle_value_type()) {
    assert(num_nodes() > 0);
    assert(has_node(a) && has_node(b) && has_node(c));
    assert(a != b && a != c && b != c);
    
    // TODO: should we sort node uids ?
    
    idx_type idx = i2u_triangles_.size();
    size_type t_uid = internal_triangles_.size();

    internal_triangles_.push_back(InternalTriangle(i2u_nodes_[a.index], i2u_nodes_[b.index], i2u_nodes_[c.index], idx, value);
    i2u_triangles_.push_back(t_uid);

    return Triangle(this, t_uid, i2u_nodes_[a.index], i2u_nodes_[b.index], i2u_nodes_[c.index]);
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
      return Triangle(m_, 
        m_->i2u_triangles_[idx_], 
        m_->internal_triangles_[idx_].uid2,
        m_->internal_triangles_[idx_].uid3 );
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
    return triangle_iterator( this, i2u_triangles_.size() );
  }







  /** @class Graph::Node2TriangleIterator
   * @brief Iterator class for triangles incident to a node. A forward iterator. */
  class Node2TriangleIterator : private totally_ordered<Node2TriangleIterator>  {

    /** Reference to the mesh */
    Mesh* m_;
    /** Node idx */
    idx_type idx_;
    /** Node idx in inner vector adj_list_n2t_ */
    idx_type idx2_;

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


    /** Construct an invalid Node2TriangleIterator. */
    Node2TriangleIterator(const Mesh* m, const Node node, idx_type idx2)
      : g_(const_cast<Mesh*>(m)) {
        assert(m_ != nullptr);
        assert(m_->num_nodes() > 0);
        assert(m_->num_triangles() > 0);
        assert(m_->has_node(node)); // the cosntructure must receive a node, instead of the idx in order to do this validation
        assert(g_->adj_list_n2t_.size() > node.index());

        idx_ = node.index();
        idx2_ = idx2;
    }

    /**
     * Reference operator for Triangle.
     * Complexity: O(1).
     *
     * @return Triangle object.
     */
    Triangle operator*() const {
      return Triangle(m_, 
        m_->i2u_triangles_[idx_], 
        m_->internal_triangles_[idx_].uid2,
        m_->internal_triangles_[idx_].uid3 );
    }

    /**
    * Node2TriangleIterator operator++. It will get to the next triangle.
    * Complexity: O(1).
    *
    * @return Node2TriangleIterator object.
    */
    Node2TriangleIterator& operator++() {
      ++idx2_;
      assert(idx2_ <= m_->adj_list_n2t_.size());
      return *this;
    }

    /**
    * Compare Node2TriangleIterator's.
    * Complexity: O(1).
    *
    * @return bool if both iterators are equal.
    */
    bool operator==(const Node2TriangleIterator& tit) const {
      return std::tie(idx_, idx2_, m_) == std::tie(tit.idx_, tit.idx2_, tit.m_);
    }

   private:
    friend class Mesh;

  };

  /**
   * Return a node2triangle_iterator pointing to the begining
   * Complexity: O(1).
   *
   * @return TriangleIterator
   */
  node2triangle_iterator node2triangle_begin(Node n) const {
    return node2triangle_iterator( this, n.index(), 0);
  }

  /**
   * Return a node2triangle_iterator pointing to one pass the last valid position.
   * Complexity: O(1).
   *
   * @return TriangleIterator
   */
  node2triangle_iterator node2triangle_end(Node n) const {
    return node2triangle_iterator( this, n.index(), adj_list_n2t_[n.index()].size() );
  }



 private:

  /** Stores all triangle objects, indexed by triangle idxs */
  std::vector<internal_triangle> internal_triangles_;

  /** Stores references to valid triangle. Vector of vector of idxs,
  indexed by node uids */
  std::vector<idx_list_type> adj_list_n2t_;

  /** Vector of adjacent triangles idxs */
  std::vector<idx_type> i2u_triangles_;

  /** Keep track of the number of triangles */
  size_type num_triangles_ = 0;

};

#endif