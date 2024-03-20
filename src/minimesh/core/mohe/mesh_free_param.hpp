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
	// parameterize by free boundary lscm
	//
	void parameterize();

	bool is_boundary(const int he_id);
	void get_pinned();
	void mass_matrix();
	void u_matrix();
	double geo_dist(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);
	void init_matrices();
	void set_matrices();
	void get_vertices(std::vector<int>& ids, std::vector<Eigen::Vector3d>& pos,int f_id, int& v1, int& v2, int& v3);
	double compute_dt(std::vector<Eigen::Vector2d>& uv);
	std::pair<bool, int> is_pinned(int& v);
	std::pair<double, double> compute_w(Eigen::Vector2d& f, Eigen::Vector2d& s);
	void A_matrix();
	void b_matrix();
	void solve_system();
	void update_positions();
	void basis_coords(std::vector<Eigen::Vector2d>& coords, const std::vector<Eigen::Vector3d>& positions);

	
private:
	// pointer to the mesh that we are working on.
	Mesh_connectivity & _m;

	// Split half edges
	std::stack<int> split_half_edges;

	// Index of boundary vertices
	// std::stack<int> boundary_ids;
	std::vector<int> boundary_ids;
	int p1;
	int p2;
    
	// map from MF1/MF2 column index to vertex id
	std::map<int, int> free;
	std::map<int, int> free_rev; //reverse

	// std::map<int, int> interior_rev;
	// map of new vertex positions
    std::map<int, Eigen::Vector3d> new_positions;

    Eigen::SparseMatrix<double> MF1;
	std::vector<Eigen::Triplet<double>> MF1_elem;
    // Eigen::SparseMatrix<double> MP1;
    Eigen::MatrixXd MP1;
	std::vector<Eigen::Triplet<double>> MP1_elem;
    Eigen::SparseMatrix<double> MF2;
	std::vector<Eigen::Triplet<double>> MF2_elem;
    // Eigen::SparseMatrix<double> MP2;
    Eigen::MatrixXd MP2;
	std::vector<Eigen::Triplet<double>> MP2_elem;

	Eigen::MatrixXd UP1;
	Eigen::MatrixXd UP2;

    Eigen::SparseMatrix<double> A;
	Eigen::MatrixXd b;
	Eigen::MatrixXd X;

};


} // end of mohe
} // end of minimesh
