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
	// std::cout << "A: " << A << std::endl;
	// std::cout << "Ubar: " << Ubar << std::endl;
	// std::cout << "Vbar: " << Vbar << std::endl;

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(A);
    Eigen::MatrixXd U = solver.solve(Ubar);
    Eigen::MatrixXd V = solver.solve(Vbar);
    // Eigen::MatrixXd U = solver.solve(Ubar).normalized();
    // Eigen::MatrixXd V = solver.solve(Vbar).normalized();

	// std::cout << "U: " << U << std::endl;
	// std::cout << "V: " << V << std::endl;

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
	double u_sum, v_sum = 0.0;
	
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


void Mesh_fixed_param::basis_coords(std::vector<Eigen::Vector2d>& coords, const std::vector<Eigen::Vector3d>& positions)
{
    if (positions.size() != 4) {
        std::cerr << "Error: Four vectors are required for computing an orthonormal basis." << std::endl;
        return;
    }

    Eigen::Vector3d x = positions[1] - positions[0];
    Eigen::Vector3d y = positions[2] - positions[0];
    Eigen::Vector3d z = positions[3] - positions[0];
    
    // Compute orthonormal basis using Gram-Schmidt process
    Eigen::Vector3d xhat = x.normalized(); // Unit vector in the direction of x
    Eigen::Vector3d zhat = (xhat.cross(y)).normalized(); // Normal to the plane formed by x and y
    Eigen::Vector3d yhat = (zhat.cross(xhat)).normalized(); // Completing the basis

    // Project positions onto the local basis and extract 2D coordinates
    for (const auto& pos : positions) {
        Eigen::Vector3d pos_local = pos - positions[0];
        double x_coord = pos_local.dot(xhat);
        double y_coord = pos_local.dot(yhat);
        coords.push_back(Eigen::Vector2d(x_coord, y_coord));
    }
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

	// // Assuming Mesh_fixed_param::basis_coords(std::vector<Eigen::Vector2d>& coords, const std::vector<Eigen::Vector3d>& positions)
	// std::vector<Eigen::Vector2d> local_basis_coords;
	// std::vector<Eigen::Vector3d> positions = {I_xyz, J_xyz, P_xyz, Q_xyz};
	// Mesh_fixed_param::basis_coords(local_basis_coords, positions);

	// Eigen::Vector2d I_xy = local_basis_coords[0];
	// Eigen::Vector2d J_xy = local_basis_coords[1];
	// Eigen::Vector2d P_xy = local_basis_coords[2];
	// Eigen::Vector2d Q_xy = local_basis_coords[3];
	// // Eigen::Vector2d I_xy = Eigen::Vector2d(I_xyz[0], I_xyz[1]);
	// // Eigen::Vector2d J_xy = Eigen::Vector2d(J_xyz[0], J_xyz[1]);
	// // Eigen::Vector2d P_xy = Eigen::Vector2d(P_xyz[0], P_xyz[1]);
	// // Eigen::Vector2d Q_xy = Eigen::Vector2d(Q_xyz[0], Q_xyz[1]);

	// Eigen::Vector2d IJ = J_xy-I_xy;
	// Eigen::Vector2d IP = P_xy-I_xy;
	// Eigen::Vector2d IQ = Q_xy-I_xy;

	Eigen::Vector3d IJ = J_xyz-I_xyz;
	Eigen::Vector3d IP = P_xyz-I_xyz;
	Eigen::Vector3d IQ = Q_xyz-I_xyz;

	// should these be 3vecs or 2vecs?
	double alpha = get_angle(IJ, IP);
	double beta = get_angle(IJ, IQ);
	double r = IJ.norm();

	// return (std::tan(alpha / 2) + std::tan(beta / 2)) / r;
	return 1.0;
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
			// printf("start\n");

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
	force_assert( mesh().check_sanity_slowly() );

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
		// if (v.is_boundary()) {
			v.data().xyz = new_positions[vid];
		// }	
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


	// while(!boundary.empty()) {
	// 	int index = boundary.top();
	// 	Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(index);
	// 	v.data().is_boundary = false;
	// 	boundary.pop();
	// }

	// assert(boundary.empty());
}

} // end of mohe
} // end of minimesh
