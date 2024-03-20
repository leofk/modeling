#include <minimesh/core/mohe/mesh_fixed_param.hpp>
#include <minimesh/core/util/assert.hpp>
#include <Eigen/Dense> 
#include <Eigen/SparseLU> 
#include <set>

#include <iostream>
#include <sstream>

namespace minimesh
{
namespace mohe
{

//
// flag all vertices that are on the boundary
//
void Mesh_fixed_param::flag_boundary()
{
	int i = 0;
	bool break_out = false;

	for(int vid = 0 ; vid < mesh().n_total_vertices() ; ++vid)
	{

		Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(vid);

		// find a vertex on the boundary
		if (ring.reset_boundary()) {
			if (break_out == false)
			{
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
					boundary_ids.push(curr_v_id);

					Mesh_connectivity::Half_edge_iterator next_he = he.next();
					curr_he_id = next_he.index();
					curr_v_id = next_he.origin().index();
					
				} while (start_v_id != curr_v_id);
				break_out = true;
			}
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


// compute position of vertex as projected onto a circle
//
void Mesh_fixed_param::generate_circle()
{

	num_boundary = static_cast<int>(boundary_ids.size());
    double angle_increment = 2 * M_PI / num_boundary;
    double angle = 0.0;

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
// compute AU = Ubar and AV = Vbar and add to new position map
//
void Mesh_fixed_param::compute_interior_pos()
{


	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    Eigen::MatrixXd U = solver.solve(Ubar);
    Eigen::MatrixXd V = solver.solve(Vbar);

    for (int i = 0; i < U.rows(); ++i) {
        new_positions[interior_rev[i]] = Eigen::Vector3d(U(i), V(i), 0.0);
    }
}

//
// compute values for A matrix per interior vertex
//
void Mesh_fixed_param::compute_A_i(int vid, int i)
{
	
	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(vid);
	do // *points to*
	{
		Mesh_connectivity::Vertex_iterator v_j = ring.half_edge().origin();
		if (!v_j.is_boundary())
		{
			int jid = v_j.index();
			int j = interior[jid];
			A_elem.push_back(Eigen::Triplet<double>(i, j, -lambda_ij(vid,jid)));
		}
	} while(ring.advance());

	A_elem.push_back(Eigen::Triplet<double>(i, i, 1.0));
}

//
// compute values for Ubar and Vbar matrix per interior vertex
//
void Mesh_fixed_param::compute_UVbar_i(int vid, int i)
{
	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(vid);
	double u_sum = 0.0;
	double v_sum = 0.0;
	
	do // *points to*
	{
		Mesh_connectivity::Vertex_iterator v_j = ring.half_edge().origin();
		if (v_j.is_boundary())
		{
			int jid = v_j.index();
			double k_ij = lambda_ij(vid,jid);

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

// get angle between 2 vectors
// double Mesh_fixed_param::get_angle(const Eigen::Vector2d &v1, const Eigen::Vector2d &v2) 
double Mesh_fixed_param::get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) 
{
	return std::acos(v1.dot(v2) / (v1.norm() * v2.norm()));
}


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
// compute mean-value weights given vertices with index i and j
//
double Mesh_fixed_param::lambda_ij(int i, int j)
{
 	Mesh_connectivity::Vertex_iterator v_i = mesh().vertex_at(i);
 	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
	double w_sum = 0.0;
	double w_ij = 0.0;

	do // *points to*
	{
		int k = ring.half_edge().origin().index();
		double w_ik = get_wik(ring.half_edge().index());

		if (k==j) { 
			w_ij = w_ik; 
		}
		
		w_sum += w_ik;

	} while(ring.advance());

	return w_ij / w_sum;
}

//
// fixed_param 
//
void Mesh_fixed_param::parametrize()
{
	flag_boundary();

	generate_circle();

	math();

	compute_interior_pos();

	update_vertex_pos();

	reset_flags();
}


//
// compute vertex updates
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

//
// reset the flags applied in subdivision step. eg. is_split
//
void Mesh_fixed_param::reset_flags() 
{
	// reset split edge flag
	while(!split_half_edges.empty()) {
		int index = split_half_edges.top();
		Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(index);
		he.data().is_split = false;
		split_half_edges.pop();
	}

	assert(split_half_edges.empty());
}

} // end of mohe
} // end of minimesh
