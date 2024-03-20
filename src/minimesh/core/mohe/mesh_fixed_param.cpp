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
				boundary_ids.push(curr_v_id);

				Mesh_connectivity::Half_edge_iterator next_he = he.next();
				curr_he_id = next_he.index();
				curr_v_id = next_he.origin().index();
				
			} while (start_v_id != curr_v_id);
			break;
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
// compute values for A matrix per interioir vertex
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
		Mesh_connectivity::Vertex_iterator v_k = ring.half_edge().origin();
		int k = v_k.index();
		double w_ik;

		Eigen::Vector3d i_pos = v_i.xyz();
		Eigen::Vector3d k_pos = v_k.xyz();
		
		Eigen::Vector3d diff = k_pos - i_pos;
   		double r_ik = diff.norm();

		double a_ik, b_ki;
		compute_angles(ring.half_edge().index(), i_pos, k_pos, a_ik, b_ki);
		w_ik = (tan(a_ik/2.0) + tan(b_ki/2.0)) / r_ik;

		w_ik = 1.0; // UNIFORM. TODO. USE MVW

		if (k==j) { 
			w_ij = w_ik; 
		}
		
		w_sum += w_ik;

	} while(ring.advance());

	return w_ij / w_sum;
}

void Mesh_fixed_param::compute_angles(int r_id, Eigen::Vector3d i_pos, Eigen::Vector3d j_pos, double& a_ik, double& b_ki)
{
	Mesh_connectivity::Half_edge_iterator r = mesh().half_edge_at(r_id); // points to i, origin at j

	Eigen::Vector3d B_pos = r.twin().next().dest().xyz();
	Eigen::Vector3d A_pos = r.next().dest().xyz();

   	Eigen::Vector3d IJ = j_pos - i_pos;
   	Eigen::Vector3d A = A_pos - i_pos;
   	Eigen::Vector3d B = B_pos - i_pos;
	double A_n = A.norm();
	double B_n = B.norm();
	double IJ_n = IJ.norm();

	a_ik = acos(A.dot(IJ) / (A_n * IJ_n));
	b_ki = acos(IJ.dot(B) / (B_n * IJ_n));
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

	int i = 0;
	// TODO - slow to iterate over all this, maybe a faster solution?
	for(int vid = 0 ; vid < mesh().n_total_vertices() ; ++vid)
	{
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(vid);

		if (!v.is_boundary())
		{
			interior[vid] = i;
			interior_rev[i] = vid;
			i++;	
		}
	}

	for (const auto& pair : interior) {
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
