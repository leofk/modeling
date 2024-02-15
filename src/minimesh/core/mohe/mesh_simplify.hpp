#pragma once

//
// mesh_simplify.hpp
//
// functionality for simplifying a mesh based on the edge-collapse schema
//
// Author: Leo FK
//


#include <string>

#include <minimesh/core/mohe/mesh_connectivity.hpp>


// Define a struct to represent an error along with its value
struct Error {
    double value; // The error value
	int v_id; // vertex associated with error

    // Constructor to initialize an Error object
    Error(double val, int v_id) : value(val), v_id(v_id) {}
    
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

class Mesh_simplify
{
public:
	// Trivial constructor
    Mesh_simplify(Mesh_connectivity& mesh_in): _m(mesh_in) {
		// printf("constructor called \n");
        // initialize_Q_matrices();
		// printf("matrices done called \n");
		// init_pos_and_errors(); 
		// printf("pos done called \n");
    }
	
	// Get the underlying mesh
	Mesh_connectivity & mesh() { return _m; }
	const Mesh_connectivity & mesh() const { return _m; }

	//
	// simplify
	//
	void simplify(const int num_entities_to_simplify);

	void initialize_Q_matrices();
	void init_pos_and_errors();
	void compute_position_and_error(Mesh_connectivity::Half_edge_iterator he);
	void mark_as_split(Mesh_connectivity::Half_edge_iterator he);
	void reset_flags();
	void collapse_edge(Mesh_connectivity::Half_edge_iterator he);
	bool check_topology(Mesh_connectivity::Half_edge_iterator he);

private:
	// pointer to the mesh that we are working on.
	Mesh_connectivity & _m;

	// A stack to split half edges
	std::stack<int> split_half_edges;

	// a map with vertex index for key and new position as value
    std::map<int, Eigen::Matrix4d> Q_matrices;

	// half edge to new vertex pos if collapsed
    std::map<int, Eigen::Vector3d> new_pos;

    std::priority_queue<Error> errorQueue;
};


} // end of mohe
} // end of minimesh
