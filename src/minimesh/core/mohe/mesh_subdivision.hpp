#pragma once

//
// mesh_subdivision.hpp
//
// functionality for subdividing a mesh based on loop's schema
//
// Author: Leo FK
//


#include <string>

#include <minimesh/core/mohe/mesh_connectivity.hpp>

namespace minimesh
{
namespace mohe
{


class Mesh_subdivision
{
public:
	// Trivial constructor
	Mesh_subdivision(Mesh_connectivity & mesh_in): _m(mesh_in) {}

	// Get the underlying mesh
	Mesh_connectivity & mesh() { return _m; }
	const Mesh_connectivity & mesh() const { return _m; }

	//
	// Given a half edge, split the edge into two
	//
	void edge_split(const int he_index);

	//
	// create additional faces. cut the first corner found, unless the face is already a triangle
	//
	void cut_a_corner(const int face_index);
	
	//
	// checks to see if a given face is a triangle
	//
	bool is_triangle(const int face_index);

	//
	// checks to see if a given half edge is at a boundary
	//
	bool is_boundary(const int he_index);

	//
	// subdivision
	//
	void subdivision();

	//
	// reset the flags applied in subdivision step. eg. is_new and is_split
	//
	void reset_flags(); 

	//
	// loop schema mask for control interior vertex
	//
	Eigen::Vector3d loop_interior_control(const int v_index); 

	//
	// loop schema mask for new interior vertex
	//
	Eigen::Vector3d loop_interior_new(const int he_index); 

	//
	// loop schema mask for control boundary vertex
	//
	Eigen::Vector3d loop_boundary_control(const int he_index);

	//
	// loop schema mask for new boundary vertex
	//
	Eigen::Vector3d loop_boundary_new(const int v_index); 


private:
	// pointer to the mesh that we are working on.
	Mesh_connectivity & _m;

	// A stack to split half edges
	std::stack<int> split_half_edges;

	// A queue to new vertices
	std::queue<int> new_vertices;

	// A queue of new precomputed positions for control/old vertices
	std::queue<Eigen::Vector3d> old_vertices_new_pos;

	// A queue of new precomputed positions for new half edge vertices
	std::queue<Eigen::Vector3d> new_vertices_pos;

	// a map with vertex index for key and new position as value
    std::map<int, Eigen::Vector3d> vertex_pos_map;

};


} // end of mohe
} // end of minimesh
