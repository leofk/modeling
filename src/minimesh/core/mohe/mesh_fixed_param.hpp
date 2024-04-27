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
	// Harmonic Parameterization
	// Fixed boundary parameterization for a open manifold mesh 
	// Fixed to the unit circle using mean-value weights
	//
	void parameterize();
	bool is_boundary(const int he_id);
	void flag_boundary();
	void generate_circle(); 
	void compute_A_i(int vid, int i); 
	void compute_UVbar_i(int vid, int i); 
	void compute_interior_pos(); 
	void update_vertex_pos(); 
	void math(); 
	double lambda_ij(const int i, const int j);
	double get_wik(int he_index);
	double get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);

private:
	// pointer to the mesh that we are working on.
	Mesh_connectivity & _m;

	// Split half edges
	std::stack<int> split_half_edges;

	// Index of boundary vertices
	std::stack<int> boundary_ids;
	int num_boundary;
    
	// vid -> mat id
	std::map<int, int> interior;

	// mat id -> vid
	std::map<int, int> interior_rev;

	// map of new vertex positions
    std::map<int, Eigen::Vector3d> new_positions;

	int RADIUS = 1;

	Eigen::MatrixXd Ubar;
	Eigen::MatrixXd Vbar;
    Eigen::SparseMatrix<double> A;
	std::vector<Eigen::Triplet<double>> A_elem;
};


} // end of mohe
} // end of minimesh
