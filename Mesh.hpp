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
  typedef typename GraphType::Node Node;
  typedef typename GraphType::Edge Edge;
  
   
  /** @struct InternalTriangle */
  struct InternalTriangle
  {
    size_type uid1;
    size_type uid2;
    size_type uid3;
    idx_type idx;
    triangle_value_type value;

    InternalTriangle(size_type uid1, size_type uid2, size_type uid3, size_type idx, triangle_value_type value) 
      : uid1(uid1), uid2(uid2), uid3(uid3), idx(idx), value(value) {
    }
  };

  /** Type of InternalTriangle */
  typedef InternalTriangle internal_triangle;

  /** Inner vector of idxs for nodes to triangles */
  //typedef std::vector<idx_type> idx_list_type;

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
    
  }


  /** Returns the number of triangles */
  size_type num_triangles() const {
    return num_triangles_();
  }


  class TNode : private totally_ordered<TNode>  {

  };

  class TEdges : private totally_ordered<TEdges>  {

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
        return Node(m_, uid1_);
      }

      /** Return node 2 of this Triangle */
      Node node2() const {
        return Node(m_, uid2_);
      }

      /** Return node 3 of this Triangle */
      Node node3() const {
        return Node(m_, uid3_);
      }



      double area() const {
        // http://www.mathopenref.com/coordtrianglearea.html
        Point a = internal_modes_[uid1_].point;
        Point b = internal_modes_[uid2_].point;
        Point c = internal_modes_[uid3_].point;
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
        return std::tie(m_, node1_uid, node2_uid, node3_uid) < std::tie(x.m_, x.node1_uid, x.node2_uid, x.node3_uid);
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
        //return  internal_triangles_[t_uid_].value;
      }

      /**
       * Returns a reference to this trinagle's value of type T as a constant.
       *
       * @return Object of type T by reference as a constant.
       */
      const triangle_value_type& value() const {
        //return  internal_triangles_[t_uid_].value;
      }


    private:
      
      Triangle(Mesh* m, size_type t_uid, size_type uid1, size_type uid2, size_type uid3) 
        : m_(const_cast<Mesh*>(m)), t_uid_(t_uid), uid1_(uid1), uid2_(uid2), uid3_(uid3) {
      }

      Mesh* m_;
      Node node1;
      Node node2;
      Node node3;
      size_type t_uid_;

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

  }





 private:

  GraphType g_nodes;
  GraphType g_triangles;

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