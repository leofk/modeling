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

class Mesh_deform
{
public:

	// constructor
    Mesh_deform(Mesh_connectivity& mesh_in): _m(mesh_in) {

	}
	
	// Get the underlying mesh
	Mesh_connectivity & mesh() { return _m; }
	const Mesh_connectivity & mesh() const { return _m; }

	//
	// Mesh Deformation
	//
	void init();
	void build_W_matrix();
	void build_L_matrix();
	void update_L_matrix();
	void deform(std::vector<int> _constraint_vids, int _handle_vid);
	void Mesh_deform::compute_r_matrices();
	void Mesh_deform::build_b_matrix();
	void Mesh_deform::solve_system();
	double Mesh_deform::compute_energy();
	void Mesh_deform::update_positions();
	double Mesh_deform::compute_wij(int he_id);
	double Mesh_deform::get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);
	bool Mesh_deform::is_constraint(int vid);

private:
	// pointer to the mesh that we are working on.
	Mesh_connectivity & _m;

    std::map<int, Eigen::MatrixXd> r_matrices;
	std::vector<int> constraint_vids;
	int handle_vid;

    Eigen::SparseMatrix<double> L;
	std::vector<Eigen::Triplet<double>> L_elem;
	Eigen::MatrixXd b;
	Eigen::MatrixXd p_prime;
	Eigen::MatrixXd p_prev;
	Eigen::MatrixXd W;

	double THRESHOLD = 1.0e-3;
};


} // end of mohe
} // end of minimesh
