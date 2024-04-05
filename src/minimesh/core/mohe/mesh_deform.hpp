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
	void clearROI() { _roi_handle.clear(); _roi_anchor.clear(); _roi_constraints.clear();};

	// append to item to region of interest
	void appendToROIHANDLE(int key, int value) {

		// reject invalid indexes
		if (key < 0) return;

		// append data to map
		_roi_handle[key] = value;
		_roi_constraints[key] = value;

		constraint_vids.push_back(value);

		// std::cout << "[" << key << "," << value << "]" << std::endl;
	}
	void appendToROIANCHOR(int key, int value) {

		// reject invalid indexes
		if (key < 0) return;

		// append data to map
		_roi_anchor[key] = value;
		_roi_constraints[key] = value;

		constraint_vids.push_back(value);

		// std::cout << "[" << key << "," << value << "]" << std::endl;
	}

	// clicked var
	int clickedVertex = -1;

	// get the region of interest
	std::map<int, int> ROI_ANCHOR() { return _roi_anchor; };
	std::map<int, int> ROI_HANDLE() { return _roi_handle; };

	//
	// Mesh Deformation
	//
	void init();
	void build_W_matrix();
	void build_L_matrix();
	void update_L_matrix();
	void deform(int _handle_id, Eigen::Vector3f pull_amount);
	void build_b_matrix();
	void build_c_map();
	void solve_system();
	double compute_energy();
	void update_positions();
	void compute_r_matrices();
	double compute_wij(int he_id);
	double get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);
	bool is_constraint(int vid);
	std::vector<int> get_constraints_vids(){ return constraint_vids; };

private:
	// pointer to the mesh that we are working on.
	Mesh_connectivity & _m;

	bool primed = false;
	std::map<int, int> _roi_handle = {
			{2,2},
	};
	std::map<int, int> _roi_anchor = {
			{5,5},
	};
	std::map<int, int> _roi_constraints = {
			{2,2},
			{5,5},
	};
	
    std::map<int, Eigen::MatrixXd> r_matrices;
	std::vector<int> constraint_vids;
	std::map<int, int> free; // vid -> matid
	std::map<int, int> free_rev; // matid -> vid
	std::map<int, int> c_map; // vid -> matid

    // Eigen::SparseMatrix<double> L;
	Eigen::MatrixXd L;

	std::vector<Eigen::Triplet<double>> L_elem;
	Eigen::MatrixXd bf;
	Eigen::MatrixXd b;
	Eigen::MatrixXd xc;
	Eigen::MatrixXd p_prime;
	Eigen::MatrixXd p_prev;
	Eigen::MatrixXd W;
	Eigen::MatrixXd Afc;
	Eigen::MatrixXd Aff;
	Eigen::Vector3d pp_handle;
	int handle_id;

	double THRESHOLD = 1.0e-3;
};


} // end of mohe
} // end of minimesh
