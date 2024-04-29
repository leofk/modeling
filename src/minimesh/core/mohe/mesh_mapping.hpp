#pragma once

#include <string>

#include <minimesh/core/mohe/mesh_connectivity.hpp>

#define M1 1
#define M2 2

struct inter_data {
    std::vector<double> bary; 
	int f_id; 
    inter_data(int f_id, std::vector<double> bary) : bary(bary), f_id(f_id) {}
};

namespace minimesh
{
namespace mohe
{

class Mesh_mapping
{
public:
	// Trivial constructor
	Mesh_mapping(Mesh_connectivity & m1, Mesh_connectivity & m2): _m1(m1), _m2(m2) {}

	// Get the underlying mesh
	Mesh_connectivity & m1() { return _m1; }
	Mesh_connectivity & m2() { return _m2; }
	Mesh_connectivity & mesh() {
		if (current_mesh == M1)
			return _m1;
		else
			return _m2;
	}
	// const Mesh_connectivity & mesh() const { return _m; }

	void build_mapping();
	void init_ISM();
	void ISM_iteration();
	inter_data trivial_map_data(int vid);
	void yep(Mesh_connectivity & mesh_in);
	inter_data find_enclosing_face(std::vector<int> faces, Eigen::Vector3d v_pos);
	Eigen::Vector3d get_pos_from_inter_data(inter_data data);
	std::vector<double> compute_bc(int f_id, Eigen::Vector3d v_pos);
	double get_wij(Eigen::Vector3d i, Eigen::Vector3d j, std::vector<Eigen::Vector3d> &verts); 


private:
	// pointer to the mesh that we are working on.
	// Mesh_connectivity & _m;
	Mesh_connectivity & _m1;
	Mesh_connectivity & _m2;

    std::map<int, int> mesh_map;
    std::map<int, inter_data> ISM_M1;
    std::map<int, inter_data> ISM_M2;
	int current_mesh;
};

double get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) {
	return std::acos(v1.dot(v2) / (v1.norm() * v2.norm()));
}


} // end of mohe
} // end of minimesh
