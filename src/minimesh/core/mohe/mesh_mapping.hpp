#pragma once

#include <string>

#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_simplify.hpp>

#define M1 1
#define M2 2

struct Inter_data {
    std::vector<double> bary; 
	int f_id; 
};

namespace minimesh
{
namespace mohe
{

class Mesh_mapping
{
public:
	Mesh_mapping(Mesh_connectivity & m1, Mesh_connectivity & m2): _m1(m1), _m2(m2) {}

	Mesh_connectivity & mesh() {
		if (current_mesh == M1)
			return _m1;
		else
			return _m2;
	}

	void switch_mesh() {
		if (current_mesh == M1)
			current_mesh = M2;
		else
			current_mesh = M1;
	}

	std::map<int, Inter_data> & current_ISM() {
		if (current_mesh == M1)
			return ISM_M1;
		else
			return ISM_M2;
	}
	std::map<int, int> & current_FTV() {
		if (current_mesh == M1)
			return m1_ftv_map;
		else
			return m2_ftv_map;
	}

	void build_mapping();
	void init_ISM();
	Inter_data ISM_iteration(HistoryEntry &entry);
	Inter_data trivial_map_data(int vid);
	Inter_data find_enclosing_face(std::vector<int> faces, Eigen::Vector3d v_pos);
	Eigen::Vector3d get_pos_from_inter_data(Inter_data data);
	std::vector<double> compute_bc(int f_id, Eigen::Vector3d v_pos);
	double get_wij(Eigen::Vector3d i, Eigen::Vector3d j, std::vector<Eigen::Vector3d> &verts); 
	void update_positions(int mesh);
	void split_vertex(Mesh_simplify& simp_m);
	double get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);

private:
	Mesh_connectivity & _m1;
	Mesh_connectivity & _m2;

	std::map<int, int> mesh_map = {
		// cow1 to horse init map
		{2348,6843}, 
		{2149,9676}, 
		{1782,9078}, 
		{1450,5922}, 
	};

    std::map<int, int> m1_ftv_map; // M1 Face to M2 Vertex Map
    std::map<int, int> m2_ftv_map; // M2 Face to M1 Vertex Map
    std::map<int, Inter_data> ISM_M1; // M1 Inter-surface Map
    std::map<int, Inter_data> ISM_M2; // M2 Inter-surface Map
	int current_mesh; // Current working mesh
};


} // end of mohe
} // end of minimesh
