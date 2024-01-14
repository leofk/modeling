#pragma once

//
// mesh_modifier.hpp
//
// Some functionality for modifying meshes.
//
// Author: Shayan Hoshyari
//


#include <string>

#include <minimesh/core/mohe/mesh_connectivity.hpp>

namespace minimesh
{
namespace mohe
{


class Mesh_modifier
{
public:
	// Trivial constructor
	Mesh_modifier(Mesh_connectivity & mesh_in): _m(mesh_in) {}

	// Get the underlying mesh
	Mesh_connectivity & mesh() { return _m; }
	const Mesh_connectivity & mesh() const { return _m; }

	//
	// Given two vertices, this function return the index of the half-edge going from v0 to v1.
	// Returns mesh::invalid_index if no half-edge exists between the two vertices.
	//
	int get_halfedge_between_vertices(const int v0, const int v1);

	//
	// Flip an edge in a mesh
	// Input: The mesh, and the index of a half-edge on the edge we wish to flip
	// Return true if operation successful, and false if operation not possible
	//
	// Assumption: mesh is all triangles
	//
	// NOTE: To see how this method works, take a look at edge-flip.svg
	//
	bool flip_edge(const int he_index);

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
	// subdivision
	//
	void subdivision();

	//
	// reset the flags applied in subdivision step. eg. is_new and is_split
	//
	void reset_flags(); 

	//
	// loop schema mask for new vertex pos for old control vertex
	//
	Eigen::Vector3d loop_mask_v(const int v_index); 

	//
	// loop schema mask for new vertex pos for split half edge
	//
	Eigen::Vector3d loop_mask_he(const int he_index); 

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

};


} // end of mohe
} // end of minimesh
