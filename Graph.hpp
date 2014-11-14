#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
#include <cassert>

#include "CS207/Util.hpp"
#include "Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 * 
 * RI(graph): All i in [0, i2u_nodes_.size) then i == internal_nodes_[i2u_nodes_[i]].idx
 * 
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {

 public:

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** Type of this graph. */
  typedef Graph graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Synonym for Point */
  typedef Point point_type;
  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;
  /** Type value for idexes */
  typedef unsigned idx_type;

  /** Type value for nodes custom data */
  typedef V node_value_type;
  /** Type value for edges custom data */
  typedef E edge_value_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;

  /** @struct InternalNode */
  struct InternalNode
  {
    point_type point;
    node_value_type value;
    idx_type idx;
    InternalNode(Point point, node_value_type value, idx_type idx) : point(point), value(value), idx(idx) {
    }
  };

  /** @struct InternalEdge */
  struct InternalEdge
  {
    size_type uid2;
    edge_value_type value;
    idx_type idx;
    InternalEdge(size_type uid2, edge_value_type value, idx_type idx) : uid2(uid2), value(value), idx(idx) {
    }
  };

  /** Type of InternalNode */
  typedef InternalNode internal_node;
  /** Type of InternalEdge */
  typedef InternalEdge internal_edge;

  /** Inner vector of internal edges for the adjacency list. */
  typedef std::vector<internal_edge> adj_list_edges;
  
  /** Inner vector of valid edges adjacency list. */
  typedef std::vector<size_type> adj_list_valid_edges;


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Construct an empty graph. */
  Graph(){
  }


  /** Default destructor */
  ~Graph() = default;

  /////////////
  // General //
  /////////////

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_nodes_.size();
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    std::cout << "Graph.clean()" << std::endl;
    internal_nodes_.clear();
    internal_edges_.clear();
    i2u_nodes_.clear();
    i2u_edges_.clear();
    num_edges_ = 0;
    assert(num_nodes() == 0);
    assert(num_edges() == 0);
  }

  /////////////////
  // GRAPH NODES //
  /////////////////

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * RI(node): uid_ == i2u_nodes_[ internal_nodes_[uid_].idx ]
   *           uid_ < internal_nodes_.size()
   *           internal_nodes_[uid_].idx < i2u_nodes_.size()
   *           uid >= 0
   *           
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node>  {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      assert(valid());
      return g_->internal_nodes_[uid_].point;
    }

    /** Return this node's position. */
    Point& position() {
      assert(valid());
      return g_->internal_nodes_[uid_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(valid());
      return g_->internal_nodes_[uid_].idx;
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& x) const {
      return std::tie(uid_, g_) == std::tie(x.uid_, x.g_);
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
      return std::tie(uid_, g_) < std::tie(x.uid_, x.g_);
    }

    /**
     * Returns a reference to this node's value of type V.
     *
     * @return Object of type V by reference
     */
    node_value_type& value() {
      assert(g_->internal_nodes_.size() > uid_);
      return g_->internal_nodes_[uid_].value;
    }

    /**
     * Returns a reference to this node's value of type V as a constant.
     *
     * @return Object of type V by reference as a constant.
     */
    const node_value_type& value() const {
      assert(g_->internal_nodes_.size() > uid_);
      return g_->internal_nodes_[uid_].value;
    }

    /**
     * Set value of type node_value_type for this node
     * @param v node_value_type
     */
    void value(node_value_type v){
      g_->internal_nodes_[uid_].value = v;
    }

    /**
     * Return the degree. The number of edges incident to this node
     * @return size_type
     */
    size_type degree() const {
      return g_->i2u_edges_[index()].size();
    }

    /**
     * Returns incident_iterator poiting to the first element.
     * @return incident_iterator
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(g_, uid_, 0);
    }

    /**
     * Returns incident_iterator poiting to one elemnt past the last valid element.
     * @return incident_iterator
     */
    incident_iterator edge_end() const {
      return IncidentIterator(g_, uid_, g_->i2u_edges_[index()].size());
    }

   private:

    // Reference to the Graph object
    Graph* g_;

    // This element's unique identification number
    size_type uid_;

    Node(const Graph* g, size_type uid)
        : g_(const_cast<Graph*>(g)), uid_(uid) {
      assert(g_ != nullptr);
    }

    bool valid() const {
      return uid_ >= 0 && uid_ < g_->internal_nodes_.size()
        && g_->internal_nodes_[uid_].idx < g_->i2u_nodes_.size()
        && g_->i2u_nodes_[g_->internal_nodes_[uid_].idx] == uid_;
    }

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
  };

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }


  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    idx_type idx = i2u_nodes_.size();
    size_type uid = internal_nodes_.size();
    
    // Create node and store it in the points vector.
    internal_nodes_.push_back( InternalNode(position, value, idx) );

    // Add the uid to the vector of valid nodes.
    i2u_nodes_.push_back(uid);

    // Push back another adjacency list for this node
    internal_edges_.push_back({});
    i2u_edges_.push_back({});

    assert(internal_edges_.size() == internal_nodes_.size());
    
    return Node(this, uid);
  }

  /** Remove a node from the graph.
   * @param[in] n Node to be removed
   * @return 1 if old has_node(n), 0 otherwise
   *
   * @post new size() == old size() - result.
   *
   * Can invalidate outstanding iterators. 
   * If old has_node(@a n), then @a n becomes invalid, as do any
   * other Node objects equal to @a n. All other Node objects remain valid.
   *
   * Complexity: O(i2u_edges_[i].size()^2)
   */
  size_type remove_node ( const Node & n){
    if( !has_node(n) )
      return 0;
    
    idx_type idx = n.index();
    size_type uid = i2u_nodes_[idx];
    
    assert( internal_nodes_[ i2u_nodes_[idx] ].idx == idx );

    // Remove all incident edges before the node is removed.
    idx_type x = 0;
    adj_list_valid_edges v = i2u_edges_[idx];
    while(x < v.size()){
      // uid of each element in v
      size_type v_uid = v[x];
    

      adj_list_valid_edges adj = i2u_edges_[ internal_nodes_[v_uid].idx ];
      
      for(size_type adj_uid : adj)
        if(adj_uid == uid){
          remove_edge(Node(this, v_uid), Node(this, adj_uid));
        }

      ++x;
    }

    // Remove the vector from i2i_edges, so i2u_nodes and i2i_edges are in sync by idx
    i2u_edges_.erase(i2u_edges_.begin() + idx);
    
    // Remove the uid from the list of valid nodes
    i2u_nodes_.erase(i2u_nodes_.begin() + idx);

    // Update the idxs for all nodes subsequent to the removed node
    x = idx;
    while(x < i2u_nodes_.size()){
      internal_nodes_[ i2u_nodes_[x] ].idx = x;
      ++x;
    }

    return 1;
  }

  /** Remove a node from the graph.
   * @param[in] n NodeIterator pointing to the node to be removed
   * @return NodeIterator
   *
   * @post new size() == old size() - result.
   *
   * Can invalidate outstanding iterators. 
   * If old has_node(@a n), then @a n becomes invalid, as do any
   * other Node objects equal to @a n. All other Node objects remain valid.
   *
   * Complexity: O(i2u_edges_[i].size()^2)
   */
  node_iterator remove_node ( node_iterator n_it ){
    remove_node(*n_it);
    return *this;
  }


  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if(this != n.g_)
      return false;

    return n.uid_ == i2u_nodes_[n.index()];
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i >= 0 && i < num_nodes());
    return Node(this, i2u_nodes_[i]);
  }

  /////////////////
  // GRAPH EDGES //
  /////////////////

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge>  {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(g_, node1_uid);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(g_, node2_uid);
    }

    /** Test whether this edge and @a x are equal.
     *
     * Equal edges are from the same graph and have the same nodes.
     */
    bool operator==(const Edge& x) const {
      return std::tie(g_, node1_uid, node2_uid) == std::tie(x.g_, x.node1_uid, x.node2_uid);
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
      return std::tie(g_, node1_uid, node2_uid) < std::tie(x.g_, x.node1_uid, x.node2_uid);
    }

    double length () const {
      return norm(node1().position() - node2().position());
    }

    /**
     * Returns a reference to this edge's value of type E.
     *
     * @return Object of type E by reference
     */
    edge_value_type& value() {
      assert(g_->internal_nodes_.size() > node1_uid);
      
      size_type min_uid = std::min( node1_uid, node2_uid );
      size_type max_uid = std::max( node1_uid, node2_uid );
      
      adj_list_edges adjl = g_->internal_edges_[min_uid];
      int z = 0;
      for(internal_edge ie : adjl){
        if(ie.uid2 == max_uid){
          return g_->internal_edges_[min_uid][z].value;
        }
        ++z;
      }

      assert(false);
    }

    /**
     * Returns a reference to this edge's value of type E as a constant.
     *
     * @return Object of type E by reference as a constant.
     */
    const edge_value_type& value() const {
      size_type min_uid = std::min( node1_uid, node2_uid );
      size_type max_uid = std::max( node1_uid, node2_uid );
      
      adj_list_edges adjl = g_->internal_edges_[min_uid];
      for(internal_edge ie : adjl)
        if(ie.uid2 == max_uid)
          return ie.value;
        
      assert(false);
    }


   private:

    Edge(const Graph* g, size_type n1, size_type n2)
        : g_(const_cast<Graph*>(g)), node1_uid(n1), node2_uid(n2) {
    }

    // Reference to the graph object
    Graph* g_;

    // Node1's uid
    size_type node1_uid;

    // Node2's uid
    size_type node2_uid;

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // Nodes a and b must be different by precondition.
    assert( !(a == b) );
    assert( a.index() != b.index() );

    size_type n1_uid = i2u_nodes_[a.index()];
    size_type n2_uid = i2u_nodes_[b.index()];

    if(has_edge(a, b))
      return Edge (this, n1_uid, n2_uid);
    
    // Insert node b uid in the node a inner vector
    internal_edges_[a.index()].push_back( internal_edge(b.index(), value, i2u_edges_[a.index()].size()) );

    // Insert node a uid in the node b inner vector
    internal_edges_[b.index()].push_back( internal_edge(a.index(), value, i2u_edges_[b.index()].size()) );

    i2u_edges_[a.index()].push_back(n2_uid);
    i2u_edges_[b.index()].push_back(n1_uid);

    assert( a.index() < i2u_nodes_.size() );
    assert( b.index() < i2u_nodes_.size() );
    
    // keep track of number of edges
    ++num_edges_;
    return Edge (this, n1_uid, n2_uid);
  }

  /** Remove an edge from the graph, or return 1 if removed, 0 otherwise.
   * Updates idx of remaining edges.
   * 
   * @pre @a a and @a b are distinct valid nodes of this graph and confirm an edge
   * @pre @a a.index() < valid_edges.size
   * @pre @a b.index() < valid_edges.size
   * 
   * @return an int, 1 if removed, 0 otherwise.
   * @post has_edge(@a a, @a b) == false
   * @post If old !has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() + 1 == old num_edges().
   *
   *
   * Complexity: O(i2u_edges_[a.index()].size() + i2u_edges_[b.index()].size())
   */
  size_type remove_edge ( const Node& a, const Node& b){
    if( !has_edge(a, b) )
      return 0;
    
    size_type n1_uid = i2u_nodes_[a.index()];
    size_type n2_uid = i2u_nodes_[b.index()];
    
    // Remove the uid from the list of valid adjacent nodes
    assert( i2u_edges_.size() > a.index() );
    assert( i2u_edges_.size() > b.index() );

    idx_type x = 0;
    while(x < i2u_edges_[a.index()].size()){
      if(i2u_edges_[a.index()][x] == n2_uid){
    
        i2u_edges_[a.index()].erase( i2u_edges_[a.index()].begin() + x);
        --num_edges_;
        break;
      }
      ++x;
    }

    x = 0;
    while(x < i2u_edges_[b.index()].size()){
      if(i2u_edges_[b.index()][x] == n1_uid){
    
        i2u_edges_[b.index()].erase( i2u_edges_[b.index()].begin() + x);
        break;
      }
      ++x;
    }

    //Update idxs
    x = a.index();
    while(x < i2u_edges_[a.index()].size()){
      internal_edges_[n1_uid][x].idx = x;
      ++x;
    }

    //Update idxs
    x = b.index();
    while(x < i2u_edges_[b.index()].size()){
      internal_edges_[n2_uid][x].idx = x;
      ++x;
    }
    return 1;
  }

  size_type remove_edge ( const Edge& e){
    return remove_edge(e.node1(), e.node2());
  }

  edge_iterator remove_edge ( edge_iterator e_it ){
    remove_edge(*e_it);
    return *this;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(i2u_edges_.size() > a.index());
    assert(i2u_edges_.size() > b.index());

    //size_type n1_uid = i2u_nodes_[a.index()];
    size_type n2_uid = i2u_nodes_[b.index()];

    adj_list_valid_edges adj_edges = i2u_edges_[a.index()];

    // Check if the edge already exists, by testing if node b is present.
    for(size_type uid2 : adj_edges)
      if( uid2 == n2_uid ){
        //std::cout << "has_edge TRUE " << a.index() << " , " << b.index() << std::endl;
        return true;
      }
    
    return false;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Make sure i is valid.
    assert(i <= num_edges_);

    // https://piazza.com/class/hyf4iomlwgj542?cid=124
    edge_iterator it = edge_begin();
    for ( ; i != 0; --i)
      ++it;

    return *it;
  }

  ///////////////
  // Iterators //
  ///////////////

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
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
      return Node(g_, g_->i2u_nodes_[idx_]);
    }

    /**
     * Incremental Operator for NodeIterator.
     * Complexity: O(1).
     *
     * @Return NodeIterator object by reference.
     */
    NodeIterator& operator++() {
      ++idx_;
      return *this;
    }

    /**
     * Equialy Operator for NodeIterator.
     * Complexity: O(1).
     *
     * @Return bool, true if both NodeIterator's are equial.
     */
    bool operator==(const NodeIterator& other) const {
      return std::tie(g_, idx_) == std::tie(other.g_, other.idx_);
    }

   private:

    NodeIterator(const Graph* g, idx_type idx)
        : g_(const_cast<Graph*>(g)), idx_(idx) {
    }

    /** Reference to the graph */
    Graph* g_;
    /** Node idx */
    idx_type idx_;

    friend class Graph;
  };

  /**
   * Return a node_iterator pointing to the begining
   * Complexity: O(1).
   *
   * @return NodeIterator
   */
  node_iterator node_begin() const {
    return NodeIterator( this, 0 );
  }

  /**
   * Return a node_iterator pointing to one pass the last valid position.
   * Complexity: O(1).
   *
   * @return NodeIterator
   */
  node_iterator node_end() const {
    return NodeIterator( this, i2u_nodes_.size() );
  }

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
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

    EdgeIterator(const Graph* g, size_type idx1_, size_type idx2_)
      : g_(const_cast<Graph*>(g)), idx1_(idx1_), idx2_(idx2_)  {
        assert(g_ != nullptr);
        assert(g_->internal_edges_.size() > 0);
        assert(g_->i2u_edges_.size() >= idx1_);
    }

    /**
     * Reference operator for Edge.
     * Complexity: O(1).
     *
     * @Return Edge object.
     */
    Edge operator*() const {
      return Edge(g_, g_->i2u_nodes_[idx1_], g_->i2u_edges_[idx1_][idx2_]);
    }

    /**
    * Incremental operator. It will get to the next edge.
    *
    * @Return edge_iterator object.
    */
    edge_iterator& operator++() {
      go_to_next();
      return *this;
    }

    /**
    * Compare Edge Iterators.
    * Complexity: O(1).
    *
    * @return bool if both iterators are equal.
    */
    bool operator==(const edge_iterator& eit) const {
      return (idx1_ == eit.idx1_ && idx2_ == eit.idx2_ ) || (idx1_ == eit.idx1_ && idx2_ == eit.idx2_ );
    }

   private:
    friend class Graph;

    

    void go_to_next() {
      
      ++idx2_;
      if(idx2_ < g_->i2u_edges_[idx1_].size())
        return;

      while(idx1_ < g_->num_nodes()){
        while(idx2_ < g_->i2u_edges_[idx1_].size()){
          if(idx2_ < g_->i2u_edges_[idx1_].size())
            return;

          ++idx2_;
        }

        // the end of the horzontal vector, move one down
        ++idx1_;
        idx2_ = 0;
      }
    }

    /** Reference to the graph */
    Graph* g_;
    /** Node1 index */
    size_type idx1_;
    /** Node2 index */
    size_type idx2_;

  };

  /**
   * Return a edge_iterator pointing to the begining
   * Complexity: O(1).
   *
   * @return EdgeIterator
   */
  edge_iterator edge_begin() const {
    return EdgeIterator( this, 0, 0);
  }

  /**
   * Return a edge_iterator pointing to one pass the last valid position.
   * Complexity: O(1).
   *
   * @return EdgeIterator
   */
  edge_iterator edge_end() const {
    return EdgeIterator ( this, num_nodes(), 0);
  }

  void echo(){
    // print internal_edges_
    std::cout << ".........................." << std::endl;
    std::cout << "internal_nodes:" << std::endl;
    int z = 0;
    for(internal_node v : internal_nodes_){
      std::cout << "uid: " << z << " -> idx: " << v.idx << " - Pos: " << v.point << std::endl;
      ++z;
    }
    std::cout << ".........................." << std::endl << std::endl;


    std::cout << "i2u_nodes:" << std::endl;
    z = 0;
    for(size_type v : i2u_nodes_){
      std::cout << "idx: " << z <<  "-> uid: " << v << std::endl;
      ++z;
    }
    std::cout << ".........................." << std::endl << std::endl;


    std::cout << "internal_edges_ " << std::endl;
    z = 0;
    for(adj_list_edges adj_edges : internal_edges_){
      std::cout << "uid: " << z << "  -> ";
      for(internal_edge v : adj_edges){
        std::cout << " [uid2: " << v.uid2 << " idx: " <<  v.idx << "] ";
      }
      ++z;
      std::cout <<std::endl ;
    }
    std::cout << ".........................." << std::endl << std::endl;


    std::cout << "i2u_edges_: " << i2u_edges_.size() << std::endl;
    z = 0;
    for(adj_list_valid_edges vv : i2u_edges_){
      std::cout << "idx: " << z << " -> ";
      for(size_type v : vv){
        std::cout << " [uid2: " << v << "]";
      }
      ++z;
      std::cout <<std::endl ;
    }
    std::cout << ".........................." << std::endl << std::endl;
  }


  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>  {

    /** Reference to the graph */
    Graph* g_;
    /** Node uid */
    size_type uid_;
    /** Inner vector iterator */
    size_type idx2_;

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


    /** Construct an invalid IncidentIterator. */
    IncidentIterator(const Graph* g, size_type uid, size_type idx2)
      : g_(const_cast<Graph*>(g)), uid_(uid), idx2_(idx2)  {
        assert(g_ != nullptr);
        assert(g_->i2u_edges_.size() > 0);
        assert(g_->internal_edges_.size() > uid_);
        assert(g_->i2u_edges_[ g_->internal_nodes_[uid_].idx ].size() >= idx2_);
    }

    /**
     * Reference operator for Edge.
     * Complexity: O(1).
     *
     * @return Edge object.
     */
    Edge operator*() const {
      // get idx for uid
      idx_type i = g_->internal_nodes_[uid_].idx;
      assert(idx2_ < g_->i2u_edges_[i].size());
      return Edge(g_, uid_, g_->i2u_edges_[i][idx2_]);
    }

    /**
    * Incremental operator. It will get to the next edge.
    *
    * @return incident_iterator object.
    */
    IncidentIterator& operator++() {
      idx_type i = g_->internal_nodes_[uid_].idx;
      adj_list_valid_edges adj_valid_edges = g_->i2u_edges_[i];
      
      size_type num_adj_nodes = adj_valid_edges.size();
      ++idx2_;

      // check if we are at the end
      if(idx2_ >= num_adj_nodes)
        idx2_ = num_adj_nodes;

      assert(idx2_ <= num_adj_nodes);
      return *this;
    }

    /**
    * Compare Incident Iterators.
    * Complexity: O(1).
    *
    * @return bool if both iterators are equal.
    */
    bool operator==(const IncidentIterator& eit) const {
      return std::tie(uid_, idx2_, g_) == std::tie(eit.uid_, eit.idx2_, eit.g_);
    }

   private:
    friend class Graph;

  };

 private:

  /** Vector with points */
  std::vector<internal_node> internal_nodes_;
  /** Adjacency List for nodes/edges */
  std::vector<adj_list_edges> internal_edges_;
  /** Store the relationship between idxs and uids */
  std::vector<size_type> i2u_nodes_;
  /** Store the adjacent valid nodes */
  std::vector<adj_list_valid_edges> i2u_edges_;
  /** Keep track of the number of edges*/
  size_type num_edges_ = 0;

};

#endif
