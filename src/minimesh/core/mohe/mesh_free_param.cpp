#include <minimesh/core/mohe/mesh_free_param.hpp>
#include <minimesh/core/util/assert.hpp>
#include <Eigen/Dense> 
#include <set>

namespace minimesh
{
namespace mohe
{

void Mesh_free_param::get_pinned()
{

	// get boundary vertex positions
	for(int he_id = 0 ; he_id < mesh().n_total_half_edges() ; ++he_id)
	{
		Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);

		if (!he.is_split()) {

			// half edge is on boundary
			if (is_boundary(he_id)) 
			{

				if (!is_boundary(he.next().index())) {
					he = he.twin();
					he_id = he.index();
				}
				int start_he_id = he_id;
				int start_v_id = he.origin().index();
				int curr_v_id = start_v_id;
				int curr_he_id = start_he_id;
				int next_v_id = curr_v_id;
				int next_he_id = curr_he_id;
				do
				{
					// boundary_ids.push(curr_v_id);

					Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(curr_he_id);
					Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(curr_v_id);
					boundary_ids.push_back(curr_v_id);
					v.data().is_boundary = true;
					Mesh_connectivity::Half_edge_iterator next_he = he.next();

					next_he_id = next_he.index();
					next_v_id = next_he.origin().index();

					curr_v_id = next_v_id;
					curr_he_id = next_he_id;
					
				} while (start_v_id != curr_v_id);
				break;
			} 
			he.data().is_split = true;
			split_half_edges.push(he.index());

			he.twin().data().is_split = true;
			split_half_edges.push(he.twin().index());
		}
	}

	// find max dist pair
    float max_distance = 0.0;
	int a, b;

    for (size_t i = 0; i < boundary_ids.size(); ++i) {
        for (size_t j = i + 1; j < boundary_ids.size(); ++j) {
			Mesh_connectivity::Vertex_iterator v_i = mesh().vertex_at(boundary_ids[i]);
			Mesh_connectivity::Vertex_iterator v_j = mesh().vertex_at(boundary_ids[j]);
			
            float d = geo_dist(v_i.xyz(), v_j.xyz());
            if (d > max_distance) {
                max_distance = d;
                pinned_ids[0] = boundary_ids[i];
				pinned_ids[1] = boundary_ids[j];
            }
        }
    }
}

void Mesh_free_param::mass_matrix()
{


	for(int f_id = 0 ; f_id < mesh().n_total_faces() ; ++f_id)
	{
		int v1, v2, v3;
		get_vertices(f_id, v1, v2, v3);
	}


}

//
// fixed_param 
//
void Mesh_free_param::parametrize()
{
	// find pinned vertices (boundary with max geodesic dist.)
	get_pinned();

	// setup M matrices ; p=2, n'=num triangles, n=num v
		// M1f - real part of free vertices: n' x n-2
		// M2f - imaginary part of free vertices: n' x n-2
		// M1p - real part of pinned vertices: n' x 2
		// M2p - imaginary part of pinned vertices: n' x 2

	mass_matrix();	
	
	// Setup U vector
		// U1p - x coord of pinned vertices: 2x1
		// U2p - y coord of pinned vertices: 2x1

	// Setup A and b
	// Solve x
	// map new positions
}

void Mesh_free_param::get_vertices(int f_id, int& v1, int& v2, int& v3)
{
	//todo
}

//
// checks to see if a given half edge is on the boundary
//
bool Mesh_free_param::is_boundary(const int he_id)
{
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);

	// see if the face for given half edge (or its twin) has negative index
	return he.twin().face().index() < 0 || he.face().index() < 0; 
}

float Mesh_free_param::geo_dist(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) 
{
    return (p1 - p2).norm();
}

void Mesh_free_param::init_matrices() 
{
	int p = 2;
	MF1.resize(mesh().n_total_faces(), mesh().n_total_vertices()-p);
	MF2.resize(mesh().n_total_faces(), mesh().n_total_vertices()-p);
	MP1.resize(mesh().n_total_faces(), p);
	MP2.resize(mesh().n_total_faces(), p);
    MF1.setZero();
    MF2.setZero();
    MP1.setZero();
    MP2.setZero();
}

} // end of mohe
} // end of minimesh
