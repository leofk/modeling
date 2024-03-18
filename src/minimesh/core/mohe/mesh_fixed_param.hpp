		#pragma once

//
// mesh_simplify.hpp
//
// functionality for fixed_paraming a mesh based on the edge-collapse schema
//
// Author: Leo FK
//


#include <string>
#include <unordered_set>
#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include "minimesh/viz/mesh_buffer.hpp"
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

namespace minimesh
{
namespace mohe
{

class Mesh_fixed_param
{
public:

	// constructor
    Mesh_fixed_param(Mesh_connectivity& mesh_in): _m(mesh_in) {}
	
	// Get the underlying mesh
	Mesh_connectivity & mesh() { return _m; }
	const Mesh_connectivity & mesh() const { return _m; }

	//
	// parametrize by fixed
	//
	void parametrize();

	bool is_boundary(const int he_id);
	
	void flag_boundary();

	void reset_flags(); 
	void generate_circle(); 
	
	void compute_A_i(int vid, int i); 
	void compute_UVbar_i(int vid, int i); 
	void compute_interior_pos(); 
	void update_vertex_pos(); 
	void math(); 
	float lambda_ij(const int i, const int j);
	void compute_angles(int r_id, Eigen::Vector3d i_pos, Eigen::Vector3d k_pos, float& a_ik, float& b_ki);

private:
	// pointer to the mesh that we are working on.
	Mesh_connectivity & _m;

	// Split half edges
	std::stack<int> split_half_edges;

	// Index of boundary vertices
	std::stack<int> boundary_ids;
    std::map<int, int> boundary;
    
	// Map between interior vertex id and position in matrices
	std::map<int, int> interior;
	std::map<int, int> interior_rev;

	// // Index of interior vertices
	// std::stack<int> interior;

	// map of new vertex positions
    std::map<int, Eigen::Vector3d> new_positions;

	int RADIUS = 1;
	// Eigen::MatrixXd A;
	Eigen::MatrixXd Ubar;
	Eigen::MatrixXd Vbar;
    Eigen::SparseMatrix<double> A;
	std::vector<Eigen::Triplet<double>> A_elem;
};


} // end of mohe
} // end of minimesh
