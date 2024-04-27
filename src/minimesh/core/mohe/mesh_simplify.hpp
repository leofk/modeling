#pragma once

//
// mesh_simplify.hpp
//
// functionality for simplifying a mesh based on the edge-collapse schema
//
// Author: Leo FK
//


#include <string>
#include <unordered_set>
#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include "minimesh/viz/mesh_buffer.hpp"
#include <minimesh/core/mohe/mesh_simplify_history.hpp>

// Define a struct to represent an error along with its value
struct Error {
    double value; // The error value
	int he_id; // half-edge associated with error

    // Constructor to initialize an Error object
    Error(double val, int he_id) : value(val), he_id(he_id) {}
    
    // Overload the comparison operator for priority queue
    bool operator<(const Error& other) const {
        // Define the comparison based on error value
        // In this example, lower error values have higher priority
        return value > other.value;
    }
};


namespace minimesh
{
namespace mohe
{

class Mesh_simplify: Mesh_simplify_history
{
public:

	// constructor
    Mesh_simplify(Mesh_connectivity& mesh_in): Mesh_simplify_history(mesh_in) {}
	
	// Get the underlying mesh
    //	Mesh_connectivity & mesh() { return _m; }
    //	const Mesh_connectivity & mesh() const { return _m; }

	// Reset method to reinitialize the state of the Mesh_simplify object
	void init() {
		initialize_Q_matrices();
		init_pos_and_errors();
	}

	//
	// simplify
	//
	void simplify(int num_entities_to_simplify);
	void simplify_to_target(int num_of_target_entities);

	//
	// compute Q matrix for each vertex in the mesh
	//
	void initialize_Q_matrices();

	//
	// initialize new position and errors for given mesh
	//
	void init_pos_and_errors();

	//
	// compute the new position and error for a half edge
	//
	void compute_position_and_error(Mesh_connectivity::Half_edge_iterator he);

	//
	// collapse an edge by creating a new vertex and deactivating/modifying the local neighbourhood
	//
	void collapse_edge(Mesh_connectivity::Half_edge_iterator he);

	// 
	// check that collapsing an edge wont result in illegal topology
	//
	bool check_connectivity(Mesh_connectivity::Half_edge_iterator he);

	//
	// check no vertices adjacent to H has 3 or less valence
	//
	bool check_valence(Mesh_connectivity::Half_edge_iterator he);
 
	// color the vertices at the top of the queue
	void color_queue(Mesh_buffer &buffer, Mesh_connectivity::Defragmentation_maps &defrag, int num_e);

	void mark_as_split(Mesh_connectivity::Half_edge_iterator he);
	void reset_flags();

    // revert simplification
    void revert_simplification();

    void construct_edge_from_history(const HistoryEntry& entry);

private:
	// pointer to the mesh that we are working on.
    //	Mesh_connectivity & _m;

	// A stack to split half edges
	std::stack<int> split_half_edges;

	// a map between vertices and their associated Q matrix
    std::map<int, Eigen::Matrix4d> Q_matrices;

	// a map between half edges and the new vector position if they were to be collapsed
    std::map<int, Eigen::Vector3d> new_pos;

	// a map to keep track of the current error for an edge
	std::map<int, double> errorMap;

	// a priority queue of half edges ordered on min error
    std::priority_queue<Error> errorQueue;
};


} // end of mohe
} // end of minimesh
