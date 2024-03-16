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

class Mesh_free_param
{
public:

	// constructor
    Mesh_free_param(Mesh_connectivity& mesh_in): _m(mesh_in) {}
	
	// Get the underlying mesh
	Mesh_connectivity & mesh() { return _m; }
	const Mesh_connectivity & mesh() const { return _m; }

	//
	// parametrize by fixed
	//
	void parametrize();

	bool is_boundary(const int he_id);
	void get_pinned();
	void mass_matrix();
	float geo_dist(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);
	void init_matrices();
	void get_vertices(int f_id, int& v1, int& v2, int& v3);
	
	
private:
	// pointer to the mesh that we are working on.
	Mesh_connectivity & _m;

	// Split half edges
	std::stack<int> split_half_edges;

	// Index of boundary vertices
	// std::stack<int> boundary_ids;
	std::vector<int> boundary_ids;
	std::vector<int> pinned_ids;
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
    Eigen::SparseMatrix<double> MF1;
    Eigen::SparseMatrix<double> MP1;
    Eigen::SparseMatrix<double> MF2;
    Eigen::SparseMatrix<double> MP2;
	std::vector<Eigen::Triplet<double>> A_elem;
};


} // end of mohe
} // end of minimesh
