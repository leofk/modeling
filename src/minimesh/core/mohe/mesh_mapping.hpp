#pragma once

#include <string>

#include <minimesh/core/mohe/mesh_connectivity.hpp>

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
	Mesh_mapping(Mesh_connectivity & mesh_in): _m(mesh_in) {}

	// Get the underlying mesh
	Mesh_connectivity & mesh() { return _m; }
	const Mesh_connectivity & mesh() const { return _m; }

	void build_mapping();
	void init_ISM();
	void ISM_iteration(int m);
	inter_data trivial_map_data(int vid);
	void yep(Mesh_connectivity & mesh_in);
	inter_data find_enclosing_face(std::vector<int> faces, Eigen::Vector3d v_pos);
	Eigen::Vector3d get_pos_from_inter_data(inter_data data);


private:
	// pointer to the mesh that we are working on.
	Mesh_connectivity & _m;
	// Mesh_connectivity & _m1;
	// Mesh_connectivity & _m2;

    std::map<int, int> mesh_map;
    std::map<int, inter_data> ISM_M1;
    std::map<int, inter_data> ISM_M2;

};


} // end of mohe
} // end of minimesh
