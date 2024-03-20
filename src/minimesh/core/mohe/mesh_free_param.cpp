#include <minimesh/core/mohe/mesh_free_param.hpp>
#include <minimesh/core/util/assert.hpp>
#include <Eigen/Dense> 
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>
#include <set>
#include <iostream>
#include <sstream>
namespace minimesh
{
namespace mohe
{

void Mesh_free_param::get_pinned()
{
	for(int vid = 0 ; vid < mesh().n_total_vertices() ; ++vid)
	{

		Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(vid);

		// find a vertex on the boundary
		if (ring.reset_boundary()) {
			Mesh_connectivity::Half_edge_iterator he = ring.half_edge();
			int he_id = he.index();

			// make sure he is on the outside of boundary for easy traversal
			if (!is_boundary(he.next().index())) 
			{
				he = he.twin();
				he_id = he.index();
			}

			int start_he_id = he_id;
			int start_v_id = he.origin().index();

			int curr_he_id = start_he_id;
			int curr_v_id = start_v_id;
			
			do
			{
				Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(curr_he_id);
				Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(curr_v_id);
				v.data().is_boundary = true;
				boundary_ids.push_back(curr_v_id);

				Mesh_connectivity::Half_edge_iterator next_he = he.next();
				curr_he_id = next_he.index();
				curr_v_id = next_he.origin().index();
				
			} while (start_v_id != curr_v_id);
			break;
		}
	}
	
	// find max dist pair
    double max_distance = 0.0;

    for (size_t i = 0; i < boundary_ids.size(); ++i) {
        for (size_t j = i + 1; j < boundary_ids.size(); ++j) {
			Mesh_connectivity::Vertex_iterator v_i = mesh().vertex_at(boundary_ids[i]);
			Mesh_connectivity::Vertex_iterator v_j = mesh().vertex_at(boundary_ids[j]);
			
            double d = geo_dist(v_i.xyz(), v_j.xyz());
            if (d > max_distance) {
                max_distance = d;
                p1 = boundary_ids[i];
				p2 = boundary_ids[j];
            }
        }
    }
}

void Mesh_free_param::basis_coords(std::vector<Eigen::Vector2d>& coords, const std::vector<Eigen::Vector3d>& positions)
{
    Eigen::Vector3d x = positions[1] - positions[0];
    Eigen::Vector3d y = positions[2] - positions[0];
    
    // compute orthonormal basis
    Eigen::Vector3d xhat = x; xhat.normalize();
    Eigen::Vector3d zhat = xhat.cross(y); 
	zhat.normalize();
    Eigen::Vector3d yhat = zhat.cross(xhat); 
	yhat.normalize();
    
    // compute coordinates in local basis
    coords.push_back(Eigen::Vector2d(0, 0));
    coords.push_back(Eigen::Vector2d(x.norm(), 0));
    coords.push_back(Eigen::Vector2d(y.dot(xhat), y.dot(yhat)));
}

void Mesh_free_param::u_matrix()
{
	UP1 = Eigen::MatrixXd::Zero(2, 1);
	UP2 = Eigen::MatrixXd::Zero(2, 1);

	Mesh_connectivity::Vertex_iterator vp1 = mesh().vertex_at(p1);
	Mesh_connectivity::Vertex_iterator vp2 = mesh().vertex_at(p2);

	std::vector<Eigen::Vector3d> pos;
	pos.push_back(vp1.xyz());
	pos.push_back(vp2.xyz());
	std::vector<Eigen::Vector2d> coords;
	basis_coords(coords, pos);

	UP1(0) = coords[0][0]; // p1.x
	UP2(0) = coords[0][1]; // p1.y
	UP1(1) = coords[1][0]; // p2.x
	UP2(1) = coords[1][1]; // p2.y
}

void Mesh_free_param::mass_matrix()
{

	init_matrices();

	int count = 0;
	int first_id = -1;

	for(int f_id = 0 ; f_id < mesh().n_total_faces() ; ++f_id)
	{
		int v1_id,  v2_id, v3_id;
		std::vector<int> ids;
		std::vector<Eigen::Vector3d> pos;
		get_vertices(ids, pos, f_id, v1_id, v2_id, v3_id);
			
		std::vector<Eigen::Vector2d> coords;
		basis_coords(coords, pos);

		// twice the area of the triangle (shouldnt be 0)
		double d_t = compute_dt(coords); 

    	for (int i = 0; i < ids.size(); ++i) {
			auto res = is_pinned(ids[i]);

			std::pair<double, double> w = {0,0};
			if (i==0) { 
				w = compute_w(coords[2], coords[1]);
			} else if (i==1) {
				w = compute_w(coords[0], coords[2]);
			} else if (i==2) {
				w = compute_w(coords[1], coords[0]);
			}
			
			if (res.first) // vid is p1 or p2
			{
				MP1_elem.push_back(Eigen::Triplet<double>(f_id, res.second, w.first/std::sqrt(d_t)));
				MP2_elem.push_back(Eigen::Triplet<double>(f_id, res.second, w.second/std::sqrt(d_t)));
			} else { 
				// check that res.second is not already given a mat_id
				
				int mat_id = free_rev[res.second];

				if (mat_id == 0 && res.second != first_id)
				{
					if (count == 0) first_id = res.second;
					mat_id = count;
					free[mat_id] = res.second;
					free_rev[res.second] = mat_id;
					count++;
				}

				MF1_elem.push_back(Eigen::Triplet<double>(f_id, mat_id, w.first/std::sqrt(d_t)));
				MF2_elem.push_back(Eigen::Triplet<double>(f_id, mat_id, w.second/std::sqrt(d_t)));
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

    std::vector<Eigen::Triplet<double>> triplets;

    // Fill the top-left corner with MF1
    for (int k = 0; k < MF1.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(MF1, k); it; ++it) {
            triplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
        }
    }

    // Fill the bottom-right corner with MF1
    for (int k = 0; k < MF1.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(MF1, k); it; ++it) {
            triplets.push_back(Eigen::Triplet<double>(it.row() + MF1.rows(), it.col() + MF1.cols(), it.value()));
        }
    }

    // Fill the top-right corner with -MF2
    for (int k = 0; k < MF2.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(MF2, k); it; ++it) {
            triplets.push_back(Eigen::Triplet<double>(it.row(), it.col() + MF1.cols(), -it.value()));
        }
    }

    // Fill the bottom-left corner with MF2
    for (int k = 0; k < MF2.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(MF2, k); it; ++it) {
            triplets.push_back(Eigen::Triplet<double>(it.row() + MF1.rows(), it.col(), it.value()));
        }
    }

    // Set the triplets into matrix A
    A.setFromTriplets(triplets.begin(), triplets.end());
}


void Mesh_free_param::b_matrix()
{
    // Construct the larger matrices by concatenating the smaller matrices
    Eigen::MatrixXd combined_MP(2 * MP1.rows(), MP1.cols() + MP2.cols());
    combined_MP << MP1, -MP2,
                   MP2, MP1;

    Eigen::MatrixXd combined_UP(UP1.rows() + UP2.rows(), UP1.cols());
    combined_UP << UP1,
                   UP2;

    // Compute the result
    Eigen::VectorXd b_temp = combined_MP * combined_UP;

    // Negate the result to obtain the final b matrix
    b = -b_temp;
}

void Mesh_free_param::solve_system()
{
	// Solve AtAx=Atb

	// initialize
    Eigen::SparseMatrix<double> At = A.transpose();
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(At * A);
    Eigen::VectorXd Atb = At * b;
    
    // backsolve
    X = solver.solve(Atb);

}

void Mesh_free_param::update_positions()
{
	// update free vertices
	// X : 2(n-p) * 1; u for all free, then v for all free
    int n = X.rows();
    int half = n / 2;

    // std::cout << "X: " << X << std::endl; 

    // Outer loop for iterating over the halves of the rows
    for (int i = 0; i < half; ++i) {
		int vid = free[i];
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(vid);
		v.data().xyz = Eigen::Vector3d(X(i), X(i+half), 0.0);
	}


	// project pinned vertices
	
	Mesh_connectivity::Vertex_iterator vp1 = mesh().vertex_at(p1);
	vp1.data().xyz = Eigen::Vector3d(UP1(0), UP2(0), 0.0);
	Mesh_connectivity::Vertex_iterator vp2 = mesh().vertex_at(p2);
	vp2.data().xyz = Eigen::Vector3d(UP1(1), UP2(1), 0.0);
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

double Mesh_free_param::compute_dt(std::vector<Eigen::Vector2d>& uv)
{
	double a = uv[0][0] * uv[1][1] - uv[0][1] * uv[1][0]; // x1y2 - y1x2
	double b = uv[1][0] * uv[2][1] - uv[1][1] * uv[2][0];  // x2y3 - y2x3
	double c = uv[2][0] * uv[0][1] - uv[2][1] * uv[0][0];  // x3y1 - y3x1
	return a+b+c;
}

std::pair<double, double> Mesh_free_param::compute_w(Eigen::Vector2d& f, Eigen::Vector2d& s)
{
	double w_real = f[0] - s[0]; // x_f - x_s
	double w_imag = f[1] - s[1]; // y_f - y_s
	return {w_real, w_imag};
}

void Mesh_free_param::get_vertices(std::vector<int>& ids,std::vector<Eigen::Vector3d>& pos, int f_id, int& v1, int& v2, int& v3)
{
	Mesh_connectivity::Face_iterator f = mesh().face_at(f_id);
	Mesh_connectivity::Half_edge_iterator he = f.half_edge();

	v1 = he.origin().index();
	Mesh_connectivity::Vertex_iterator vv1 = mesh().vertex_at(v1);
	v2 = he.dest().index();
	Mesh_connectivity::Vertex_iterator vv2 = mesh().vertex_at(v2);
	v3 = he.next().dest().index();
	Mesh_connectivity::Vertex_iterator vv3 = mesh().vertex_at(v3);


	ids.push_back(v1);
	ids.push_back(v2);
	ids.push_back(v3);

	pos.push_back(vv1.xyz());
	pos.push_back(vv2.xyz());
	pos.push_back(vv3.xyz());
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

double Mesh_free_param::geo_dist(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) 
{
    return (p1 - p2).norm();
}

void Mesh_free_param::init_matrices() 
{
	int p = 2;
	MF1.resize(mesh().n_total_faces(), mesh().n_total_vertices()-p);
	MF2.resize(mesh().n_total_faces(), mesh().n_total_vertices()-p);
	// MP1.resize(mesh().n_total_faces(), p);
	// MP2.resize(mesh().n_total_faces(), p);
    MF1.setZero();
    MF2.setZero();
    MP1.setZero();
    MP2.setZero();

	MP1 = Eigen::MatrixXd::Zero(mesh().n_total_faces(), 2);
    MP2 = Eigen::MatrixXd::Zero(mesh().n_total_faces(), 2);


	// A.resize(2 * MF1.rows(), 2 * MF1.cols());
    // A.setZero();
}

void Mesh_free_param::set_matrices() 
{
	MF1.setFromTriplets(MF1_elem.begin(), MF1_elem.end());
	MF2.setFromTriplets(MF2_elem.begin(), MF2_elem.end());
	// MP1.setFromTriplets(MP1_elem.begin(), MP1_elem.end());
	// MP2.setFromTriplets(MP2_elem.begin(), MP2_elem.end());

    for (const auto& triplet : MP1_elem) {
        MP1(triplet.row(), triplet.col()) = triplet.value();
    }

    for (const auto& triplet : MP2_elem) {
        MP2(triplet.row(), triplet.col()) = triplet.value();
    }
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
