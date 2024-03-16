#include <minimesh/core/mohe/mesh_free_param.hpp>
#include <minimesh/core/util/assert.hpp>
#include <Eigen/Dense> 
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>
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

    for (size_t i = 0; i < boundary_ids.size(); ++i) {
        for (size_t j = i + 1; j < boundary_ids.size(); ++j) {
			Mesh_connectivity::Vertex_iterator v_i = mesh().vertex_at(boundary_ids[i]);
			Mesh_connectivity::Vertex_iterator v_j = mesh().vertex_at(boundary_ids[j]);
			
            float d = geo_dist(v_i.xyz(), v_j.xyz());
            if (d > max_distance) {
                max_distance = d;
                p1 = boundary_ids[i];
				p2 = boundary_ids[j];
            }
        }
    }
}

void Mesh_free_param::u_matrix()
{
	UP1 = Eigen::MatrixXd::Zero(2, 1);
	UP2 = Eigen::MatrixXd::Zero(2, 1);

	Mesh_connectivity::Vertex_iterator vp1 = mesh().vertex_at(p1);
	Mesh_connectivity::Vertex_iterator vp2 = mesh().vertex_at(p2);

	UP1(0) = vp1.xyz()[0]; // p1.x
	UP2(0) = vp1.xyz()[1]; // p1.y
	UP1(1) = vp2.xyz()[0]; // p2.x
	UP2(1) = vp2.xyz()[1]; // p2.y
}

void Mesh_free_param::mass_matrix()
{

	init_matrices();

	int i = 0;

	for(int f_id = 0 ; f_id < mesh().n_total_faces() ; ++f_id)
	{
		int v1_id, v2_id, v3_id;
		std::vector<int> ids;
		get_vertices(ids, f_id, v1_id, v2_id, v3_id);

		float d_t = compute_dt(v2_id, v2_id, v3_id); 
		auto w = compute_w(v2_id, v3_id);

    	for (int i = 0; i < ids.size(); ++i) {
			auto res = is_pinned(ids[i]);
			if (res.first) // vid is p1 or p2
			{
				MP1_elem.push_back(Eigen::Triplet<double>(f_id, res.second, w.first/d_t));
				MP2_elem.push_back(Eigen::Triplet<double>(f_id, res.second, w.second/d_t));
			} else { 
				free[i] = res.second;
				MF1_elem.push_back(Eigen::Triplet<double>(f_id, i, w.first/d_t));
				MF2_elem.push_back(Eigen::Triplet<double>(f_id, i, w.second/d_t));
				i++;
			}
		}
	}

	set_matrices();
}

void Mesh_free_param::A_matrix()
{
    // Construct matrix A
    A.resize(2 * MF1.rows(), 2 * MF1.cols());
    A.reserve(2 * MF1.nonZeros() + 2 * MF2.nonZeros());
    A.topLeftCorner(MF1.rows(), MF1.cols()) = MF1;
    A.bottomRightCorner(MF1.rows(), MF1.cols()) = MF1;
    A.topRightCorner(MF2.rows(), MF2.cols()) = -MF2;
    A.bottomLeftCorner(MF2.rows(), MF2.cols()) = MF2;
}

void Mesh_free_param::b_matrix()
{
    b = -(Eigen::kroneckerProduct(MP1, Eigen::MatrixXd::Identity(UP1.cols(), UP1.cols())) * UP1 +
        Eigen::kroneckerProduct(MP2, Eigen::MatrixXd::Identity(UP2.cols(), UP2.cols())) * UP2);
}

void Mesh_free_param::solve_system()
{
	// Solve ax=b
  	Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
	solver.compute(A);
	X = solver.solve(b);
}

void Mesh_free_param::update_positions()
{
	// update free vertices
	// X : 2(n-p) * 1; u for all free, then v for all free
    int num_rows = X.rows();
    int half_rows = num_rows / 2;

    // Outer loop for iterating over the halves of the rows
    for (int half = 0; half < 2; ++half) {

        // Inner loop for iterating over the rows within each half
        for (int i = half * half_rows; i < (half + 1) * half_rows; ++i) {
			int vid = free[half];
			Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(vid);
			v.data().xyz = Eigen::Vector3d(X(half), X(i), 0.0);
        }
    }

	// project pinned vertices
	Mesh_connectivity::Vertex_iterator vp1 = mesh().vertex_at(p1);
	vp1.data().xyz = Eigen::Vector3d(vp1.xyz()[0], vp1.xyz()[1], 0.0);
	Mesh_connectivity::Vertex_iterator vp2 = mesh().vertex_at(p2);
	vp2.data().xyz = Eigen::Vector3d(vp2.xyz()[0], vp2.xyz()[1], 0.0);
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
	u_matrix();

	// Setup A and b
	A_matrix();
	b_matrix();

	// Solve x
	solve_system();

	// map new positions
	update_positions();
}

float Mesh_free_param::compute_dt(int& v1_id, int& v2_id, int& v3_id)
{
	Mesh_connectivity::Vertex_iterator v1 = mesh().vertex_at(v1_id);
	Mesh_connectivity::Vertex_iterator v2 = mesh().vertex_at(v2_id);
	Mesh_connectivity::Vertex_iterator v3 = mesh().vertex_at(v3_id);

	float a = v1.xyz()[0] * v2.xyz()[1] - v1.xyz()[1] * v2.xyz()[0]; // x1y2 - y1x2
	float b = v2.xyz()[0] * v3.xyz()[1] - v2.xyz()[1] * v3.xyz()[0]; // x2y3 - y2x3
	float c = v3.xyz()[0] * v1.xyz()[1] - v3.xyz()[1] * v1.xyz()[0]; // x3y1 - y3x1

	return a+b+c;
}

std::pair<float, float> Mesh_free_param::compute_w(int& v2_id, int& v3_id)
{
	Mesh_connectivity::Vertex_iterator v2 = mesh().vertex_at(v2_id);
	Mesh_connectivity::Vertex_iterator v3 = mesh().vertex_at(v3_id);
	float w_real = v3.xyz()[0] - v2.xyz()[0]; // x_3 - x_2
	float w_imag = v3.xyz()[1] - v2.xyz()[1]; // y_3 - y_2
	return {w_real, w_imag};
}

void Mesh_free_param::get_vertices(std::vector<int>& ids, int f_id, int& v1, int& v2, int& v3)
{
	Mesh_connectivity::Face_iterator f = mesh().face_at(f_id);
	Mesh_connectivity::Half_edge_iterator he = f.half_edge();

	v1 = he.origin().index();
	v2 = he.dest().index();
	v3 = he.next().dest().index();

	ids.push_back(v1);
	ids.push_back(v2);
	ids.push_back(v3);
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

void Mesh_free_param::set_matrices() 
{
	MF1.setFromTriplets(MF1_elem.begin(), MF1_elem.end());
	MF2.setFromTriplets(MF2_elem.begin(), MF2_elem.end());
	MP1.setFromTriplets(MP1_elem.begin(), MP1_elem.end());
	MP2.setFromTriplets(MP2_elem.begin(), MP2_elem.end());
}

std::pair<bool, int> Mesh_free_param::is_pinned(int& v) {

	if (v == p1) {
		return {true, 0};
	}
	if (v == p2) {
		return {true, 1};
	}
	return {false, v};
}

} // end of mohe
} // end of minimesh
