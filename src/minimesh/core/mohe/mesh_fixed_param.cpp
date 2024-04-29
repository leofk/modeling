#include <minimesh/core/mohe/mesh_fixed_param.hpp>
#include <minimesh/core/util/assert.hpp>
#include <Eigen/Dense> 
#include <Eigen/SparseLU> 
#include <set>

#include <iostream>
#include <sstream>
#include "minimesh/core/util/numbers.hpp"

namespace minimesh
{
namespace mohe
{

//
// identify all the boundary and interior vertices on the mesh
//
void Mesh_fixed_param::flag_boundary()
{
	int i = 0; 
	bool break_out = false;

	for(int vid = 0 ; vid < mesh().n_total_vertices() ; ++vid)
	{
		Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(vid);

		// find the first vertex that lies on the boundary
		if (ring.reset_boundary()) {

			// check to see if we have already traversed the boundary
			if (break_out == false)
			{
				Mesh_connectivity::Half_edge_iterator he = ring.half_edge();
				int he_id = he.index();

				// make sure HE is on the outside of boundary for easy traversal
				if (!is_boundary(he.next().index())) 
				{
					he = he.twin();
					he_id = he.index();
				}

				int start_he_id = he_id;
				int start_v_id = he.origin().index();

				int curr_he_id = start_he_id;
				int curr_v_id = start_v_id;
				
				// traverse the boundary and tabulate an ordered list of boundary vertex data
				do
				{
					Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(curr_he_id);
					Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(curr_v_id);
					v.data().is_boundary = true;
					boundary_ids.push(curr_v_id);

					Mesh_connectivity::Half_edge_iterator next_he = he.next();
					curr_he_id = next_he.index();
					curr_v_id = next_he.origin().index();
					
				} while (start_v_id != curr_v_id);
				break_out = true;
			}

		// vertex is not on boundary, tabulate as interior
		} else {		
			interior[vid] = i;
			interior_rev[i] = vid;
			i++;	
		}
	}
}

//
// checks to see if a given half edge is on the boundary
//
bool Mesh_fixed_param::is_boundary(const int he_id)
{
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);

	// see if the face for given half edge (or its twin) has negative index
	return he.twin().face().index() < 0 || he.face().index() < 0; 
}

//
// compute positions for boundary vertices on the unit circle
//
void Mesh_fixed_param::generate_circle()
{

	num_boundary = static_cast<int>(boundary_ids.size());
    double angle_increment = 2 * numbers::pi / num_boundary;
    double angle = 0.0;

	// iterate over ordered stack of boundary vertices and compute their position on the unit circle
    while (!boundary_ids.empty()) {
        int vid = boundary_ids.top();

        double x_new = RADIUS * cos(angle);
        double y_new = RADIUS * sin(angle);

        new_positions[vid] = Eigen::Vector3d(x_new, y_new, 0.0);

        angle += angle_increment;
		
		boundary_ids.pop();
	}
}

//
// compute UVs by solving linear system Au = Ubar & Av = Vbar
//
void Mesh_fixed_param::compute_interior_pos()
{
	// Use LU since A is not symmetric!
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    Eigen::MatrixXd U = solver.solve(Ubar);
    Eigen::MatrixXd V = solver.solve(Vbar);

	// populate vertex position map with solved UVs
    for (int i = 0; i < U.rows(); ++i) {
        new_positions[interior_rev[i]] = Eigen::Vector3d(U(i), V(i), 0.0);
    }
}

//
// populate A matrix
// vid: vertex index in mesh data structure
// i: vertex index in interior matrices (A, U, V, etc.) 
//
void Mesh_fixed_param::compute_A_i(int vid, int i)
{
	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(vid);
	do // *points to*
	{
		Mesh_connectivity::Vertex_iterator v_j = ring.half_edge().origin();

		// check adjacent vertex is NOT on the boundary
		if (!v_j.is_boundary())
		{
			int jid = v_j.index();
			int j = interior[jid];
			
			// A_{ij} is inverse normalized weight for ij
			A_elem.push_back(Eigen::Triplet<double>(i, j, -lambda_ij(vid,jid)));
		}
	} while(ring.advance());

	// A_{ii} is 1.0
	A_elem.push_back(Eigen::Triplet<double>(i, i, 1.0));
}

//
// populate Ubar and Vbar vectors
// vid: vertex index in mesh data structure
// i: vertex index in interior matrices (A, U, V, etc.) 
//
void Mesh_fixed_param::compute_UVbar_i(int vid, int i)
{
	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(vid);
	double u_sum = 0.0;
	double v_sum = 0.0;
	
	do // *points to*
	{
		Mesh_connectivity::Vertex_iterator v_j = ring.half_edge().origin();
		
		// check adjacent vertex IS on the boundary
		if (v_j.is_boundary())
		{
			int jid = v_j.index();

			// normalized weight for ij
			double k_ij = lambda_ij(vid,jid);

			// new unit circle positions for j
			double u_j = new_positions[jid][0];
			double v_j = new_positions[jid][1];

			double u_ij = k_ij * u_j;
			double v_ij = k_ij * v_j;

			u_sum += u_ij;
			v_sum += v_ij;	
		}
	} while(ring.advance());

	Ubar(i) = u_sum;
	Vbar(i) = v_sum;
}

// 
// compute angle between two vectors (cosine similarity)
//
double Mesh_fixed_param::get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) 
{
	return std::acos(v1.dot(v2) / (v1.norm() * v2.norm()));
}

//
// compute mean value weights for a half edge
//
double Mesh_fixed_param::get_wik(int he_index) 
{
	// points to i from j
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_index); 

	Mesh_connectivity::Vertex_iterator I = he.dest();
	Mesh_connectivity::Vertex_iterator J = he.origin();  
	Mesh_connectivity::Vertex_iterator P = he.twin().next().dest(); 
	Mesh_connectivity::Vertex_iterator Q = he.next().dest(); 

	Eigen::Vector3d I_xyz = I.xyz();
	Eigen::Vector3d J_xyz = J.xyz();
	Eigen::Vector3d P_xyz = P.xyz();
	Eigen::Vector3d Q_xyz = Q.xyz();

	Eigen::Vector3d IJ = J_xyz-I_xyz;
	Eigen::Vector3d IP = P_xyz-I_xyz;
	Eigen::Vector3d IQ = Q_xyz-I_xyz;

	double alpha = get_angle(IJ, IP);
	double beta = get_angle(IJ, IQ);
	double r = IJ.norm();

	return (std::tan(alpha / 2) + std::tan(beta / 2)) / r;
}

//
// compute lambda_{ij} by normalizing mean-value weights between all vertices adjacent to i
//
double Mesh_fixed_param::lambda_ij(int i, int j)
{
 	Mesh_connectivity::Vertex_iterator v_i = mesh().vertex_at(i);
 	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
	double w_sum = 0.0;
	double w_ij = 0.0;

	// compute MVW for all adjacent vertices
	do // *points to*
	{
		int k = ring.half_edge().origin().index();
		double w_ik = get_wik(ring.half_edge().index());

		// k is target adajcent vertex j
		if (k==j) { 
			w_ij = w_ik; 
		}
		
		w_sum += w_ik;

	} while(ring.advance());

	return w_ij / w_sum;
}

//
// Harmonic Parameterization
// Fixed boundary parameterization for a open manifold mesh 
// Fixed to the unit circle using mean-value weights
//
void Mesh_fixed_param::parameterize()
{
	flag_boundary();

	generate_circle();

	math();

	compute_interior_pos();

	update_vertex_pos();
}


//
// setup matrices for linear system
//
void Mesh_fixed_param::math() 
{	
	int n = mesh().n_total_vertices() - num_boundary;
	A.resize(n, n);
    A.setZero();
	Ubar = Eigen::MatrixXd::Zero(n, 1);
	Vbar = Eigen::MatrixXd::Zero(n, 1);
	for (const auto& pair : interior) 
	{
		compute_A_i(pair.first, pair.second);
		compute_UVbar_i(pair.first, pair.second);
    }
	A.setFromTriplets(A_elem.begin(), A_elem.end());
}

//
// update all vertex positions
//
void Mesh_fixed_param::update_vertex_pos() 
{
	for(int vid = 0 ; vid < mesh().n_total_vertices() ; ++vid)
	{
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(vid);
		v.data().xyz = new_positions[vid];
	}
}

} // end of mohe
} // end of minimesh
