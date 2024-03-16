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
					boundary_ids.push(curr_v_id);

					Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(curr_he_id);
					Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(curr_v_id);
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

	int size = boundary_ids.size();

    double angle_increment = 2 * M_PI / size;
    double angle = 0.0;

    while (!boundary_ids.empty()) {
        int vid = boundary_ids.top();

        // Compute new xy position of vertex on the circle using trigonometry
        float x_new = RADIUS * cos(angle);
        float y_new = RADIUS * sin(angle);
		// printf("x new = %f , y new = %f \n", x_new, y_new);

        new_positions[vid] = Eigen::Vector3d(x_new, y_new, 0.0);

        // Increase the angle for the next vertex
        angle += angle_increment;
		
		boundary_ids.pop();
	}
}



//
// compute position of vertex as projected onto a circle
//
void Mesh_fixed_param::compute_circle_pos(const int vid)
{
	Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(vid);
	
	float x = v.xyz()[0];
	float y = v.xyz()[1];
	// printf("x = %f , y = %f \n", x, y);


    float theta = std::atan2(y, x); 

	if (theta < 0) {
		theta += 2 * M_PI;
	}

    float x_new = RADIUS * std::cos(theta); 
    float y_new = RADIUS * std::sin(theta);
	// printf("x new = %f , y new = %f \n", x_new, y_new);

	new_positions[vid] = Eigen::Vector3d(x_new, y_new, 0.0);
}

//
// compute AU = Ubar and AV = Vbar and add to new position map
//
void Mesh_fixed_param::compute_interior_pos()
{
    // Eigen::MatrixXd U = A.partialPivLu().solve(Ubar);
    // Eigen::MatrixXd V = A.partialPivLu().solve(Vbar);

	std::cout << "A: " << A << std::endl;
	std::cout << "Ubar: " << Ubar << std::endl;
	std::cout << "Vbar: " << Vbar << std::endl;

  	Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
	solver.compute(A);
	// TODO WORKS IF NORMALIZED??A
	// Eigen::MatrixXd U = solver.solve(Ubar).normalized();
	// Eigen::MatrixXd V = solver.solve(Vbar).normalized();
	Eigen::MatrixXd U = solver.solve(Ubar);
	Eigen::MatrixXd V = solver.solve(Vbar);


	std::cout << "U: " << U << std::endl;
	std::cout << "V: " << V << std::endl;

    for (int i = 0; i < U.rows(); ++i) {
		float u = U(i);
		float v = V(i);
		// printf("u = %f , v = %f \n", u, v);
        new_positions[interior_rev[i]] = Eigen::Vector3d(U(i), V(i), 0.0);
    }
}

//
// compute values for A matrix per interioir vertex
//
void Mesh_fixed_param::compute_A_i(int vid, int i)
{
	// A_elem.push_back(Eigen::Triplet<double>(i, i, 1.0));
	float sum = 0.0;
	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(vid);
	do // *points to*
	{
		Mesh_connectivity::Vertex_iterator v_j = ring.half_edge().origin();
		if (!v_j.is_boundary())
		{
			int j = interior[v_j.index()];
			// A(i, j) = -lambda_ij(vid,v_j.index());
			A_elem.push_back(Eigen::Triplet<double>(i, j, -lambda_ij(vid,v_j.index())));
			sum += 1.0;
			//sum += lambda_ij(vid,v_j.index());
			// printf("aij = %f \n", A(i,j));

		}
	} while(ring.advance());

	A_elem.push_back(Eigen::Triplet<double>(i, i, sum));


}

//
// compute values for Ubar and Vbar matrix per interioir vertex
//
void Mesh_fixed_param::compute_UVbar_i(int vid, int i)
{
	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(vid);
	float u_sum, v_sum = 0.0;
	
	do // *points to*
	{
		Mesh_connectivity::Vertex_iterator v_j = ring.half_edge().origin();
		if (v_j.is_boundary())
		{
			int j = interior[v_j.index()];
			u_sum += lambda_ij(vid,v_j.index()) * new_positions[v_j.index()][0]; // k_ij * u_j
			v_sum += lambda_ij(vid,v_j.index()) * new_positions[v_j.index()][1];
		}
	} while(ring.advance());
	// printf("ui = %f \n", Ubar(i));
	// printf("vi = %f \n", Vbar(i));

	Ubar(i) = u_sum;
	Vbar(i) = v_sum;
}

//
// compute mean-value weights given vertices with index i and j
//
float Mesh_fixed_param::lambda_ij(int i, int j)
{
 	Mesh_connectivity::Vertex_iterator v_i = mesh().vertex_at(i);
 	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
	float w_sum = 0.0;
	float w_ij = 0.0;

	do // *points to*
	{
		Mesh_connectivity::Vertex_iterator v_k = ring.half_edge().origin();
		int k = v_k.index();
		float w_ik;

		Eigen::Vector3d i_pos = v_i.xyz();
		Eigen::Vector3d k_pos = v_k.xyz();
		
		Eigen::Vector3d diff = k_pos - i_pos;
   		float r_ik = diff.norm();

		float a_ik, b_ki;
		compute_angles(ring.half_edge().index(), i_pos, k_pos, a_ik, b_ki);
		w_ik = (tan(a_ik/2.0) + tan(b_ki/2.0)) / r_ik;

		if (k==j) { w_ij = w_ik; }
		
		w_sum += w_ik;

	} while(ring.advance());

	//return w_ij / w_sum;
	return 1.0;
}

void Mesh_fixed_param::compute_angles(int r_id, Eigen::Vector3d i_pos, Eigen::Vector3d k_pos, float& a_ik, float& b_ki)
{
	Mesh_connectivity::Half_edge_iterator r = mesh().half_edge_at(r_id); // points to i, origin at k

	Eigen::Vector3d A_pos = r.twin().next().dest().xyz();
	Eigen::Vector3d B_pos = r.next().dest().xyz();

   	Eigen::Vector3d IJ = k_pos - i_pos;
   	Eigen::Vector3d A = A_pos - i_pos;
   	Eigen::Vector3d B = B_pos - i_pos;
	float A_n = A.norm();
	float B_n = B.norm();
	float IJ_n = IJ.norm();

	a_ik = acos(A.dot(IJ) / (A_n * IJ_n));
	b_ki = acos(B.dot(IJ) / (B_n * IJ_n));
}


//
// fixed_param 
//
void Mesh_fixed_param::parametrize()
{
	// iterated over vertices to find boundary/nonboundary vertices
	flag_boundary();
	printf("flags done \n");
	force_assert( mesh().check_sanity_slowly() );
	
	math();
	printf("math done \n");

	force_assert( mesh().check_sanity_slowly() );


	compute_interior_pos();
	printf("interior done \n");

	force_assert( mesh().check_sanity_slowly() );


	update_vertex_pos();
	printf("update done \n");

	force_assert( mesh().check_sanity_slowly() );


	reset_flags();
}


//
// compute vertex updates
//
void Mesh_fixed_param::math() 
{
	printf("start");
	
	// int n = mesh().n_total_vertices() - boundary.size();
	int n = mesh().n_total_vertices() - boundary_ids.size();
	// A = Eigen::MatrixXd::Zero(n, n);
	A.resize(n, n);
    A.setZero();
	Ubar = Eigen::MatrixXd::Zero(n, 1);
	Vbar = Eigen::MatrixXd::Zero(n, 1);

	printf("y%d", n);

    // for (const auto& pair : boundary) {
	// 	compute_circle_pos(pair.second);
    // }

	generate_circle();

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
