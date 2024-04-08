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


	// clear the region of interest
	void clear_constraints() { _handle.clear(); _fixed.clear(); _constraints.clear();};

	// append to item to region of interest
	void append_handle(int key, int value) {

		// reject invalid indexes
		if (key < 0) return;

		// append data to map
		_handle[key] = value;
		_constraints[key] = value;
		// handle_id = value;

		// constraint_vids.push_back(value);

		// std::cout << "[" << key << "," << value << "]" << std::endl;
	}

	void append_fixed(int key, int value) {

		// reject invalid indexes
		if (key < 0) return;

		// append data to map
		_fixed[key] = value;
		_constraints[key] = value;

		// constraint_vids.push_back(value);

		// std::cout << "[" << key << "," << value << "]" << std::endl;
	}

	bool is_handle(int id) { 
		return _handle.count(id) > 0;
	};

	void needs_reinit() {
		primed = false;
	}
	// clicked var
	int clickedVertex = -1;

	// get the region of interest
	std::map<int, int> fixed_map() { return _fixed; };
	std::map<int, int> handle_map() { return _handle; };

	//
	// Mesh Deformation
	//
	void init();
	void build_W_matrix();
	void build_L_matrix(Eigen::SparseMatrix<double> &Aff);
	void update_L_matrix();
	void deform(int _handle_id, Eigen::Vector3f pull_amount);
	void build_b_matrix();
	void build_c_map();
	void solve_system();
	double compute_energy();
	void update_positions();
	void compute_r_matrices();
	void compute_r_i(int vid);
	double compute_wij(int he_id);
	double get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);
	bool is_constraint(int vid);
	bool is_boundary_edge(int he_index);

private:
	// pointer to the mesh that we are working on.
	Mesh_connectivity & _m;

	bool primed = false;
	bool xc_done = false;
	bool first = true;

	std::map<int, int> _handle = {
			// {22,22}, 
			// {0,0}, 
			// {11,11}, 
	};
	std::map<int, int> _fixed = {
			// {35,35}, 
			// {46,46},
			// {0,0}, 
	};
	std::map<int, int> _constraints = {
			// {35,35}, 
			// {46,46}, 
			// {22,22}, 
			// {0,0}, 
			// {11,11}, 
	};
	// int handle_id = 22;
	
    std::map<int, Eigen::MatrixXd> r_matrices;
    std::map<int, Eigen::MatrixXd> I_matrices;
	// std::vector<int> constraint_vids;
	std::map<int, int> free; // vid -> matid
	std::map<int, int> free_rev; // matid -> vid
	std::map<int, int> c_map; // vid -> matid

    // Eigen::SparseMatrix<double> Aff;
	// Eigen::MatrixXd L;
	// std::vector<Eigen::Triplet<double>> Aff_elem;

	Eigen::MatrixXd bf;
	// Eigen::MatrixXd b;
	Eigen::MatrixXd xc;
	Eigen::MatrixXd p_prime;
	// Eigen::MatrixXd p_prev;
	Eigen::MatrixXd W;
	Eigen::MatrixXd Afc;
	// Eigen::MatrixXd Aff;
	Eigen::Vector3d pp_handle;
	std::map<int, Eigen::Vector3d> pp_handles;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

	double THRESHOLD = 1.0e-3;
};


} // end of mohe
} // end of minimesh
