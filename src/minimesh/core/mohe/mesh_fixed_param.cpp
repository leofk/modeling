#include <minimesh/core/mohe/mesh_fixed_param.hpp>
#include <minimesh/core/util/assert.hpp>
#include <Eigen/Dense> 
#include <set>

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
				he.origin().data().is_boundary = true;
				// boundary.push(he.origin().index());
				boundary[he.origin().index()] = he.origin().index();

				he.dest().data().is_boundary = true;
				// boundary.push(he.dest().index());
				boundary[he.dest().index()] = he.dest().index();
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

//
// compute position of vertex as projected onto a circle
//
void Mesh_fixed_param::compute_circle_pos(const int vid)
{
	Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(vid);
	
	float x = v.xyz()[0];
	float y = v.xyz()[1];
	// printf("x=%f , y=%f \n", x, y);


    float theta = std::atan2(y, x); 

    float x_new = RADIUS * std::cos(theta); 
    float y_new = RADIUS * std::sin(theta);
	// TODO -  maybe need to adjust based on quadrant??

	new_positions[vid] = Eigen::Vector3d(x_new, y_new, 0.0);
}

//
// compute AU = Ubar and AV = Vbar and add to new position map
//
void Mesh_fixed_param::compute_interior_pos()
{
    Eigen::MatrixXd U = A.partialPivLu().solve(Ubar);
    Eigen::MatrixXd V = A.partialPivLu().solve(Vbar);

    for (int i = 0; i < U.rows(); ++i) {
        new_positions[i] = Eigen::Vector3d(U(i), V(i), 0.0);
    }
}

//
// compute values for A matrix per interioir vertex
//
void Mesh_fixed_param::compute_A_i(const int i)
{
	A(i, i) = 1.0;

	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
	do // *points to*
	{
		Mesh_connectivity::Vertex_iterator v_j = ring.half_edge().origin();
		if (!v_j.is_boundary())
		{
			int j = v_j.index();
			A(i, j) = -lambda_ij(i,j);
		}
	} while(ring.advance());

}

//
// compute values for Ubar and Vbar matrix per interioir vertex
//
void Mesh_fixed_param::compute_UVbar_i(const int i)
{
	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
	float u_sum, v_sum = 0.0;
	
	do // *points to*
	{
		Mesh_connectivity::Vertex_iterator v_j = ring.half_edge().origin();
		if (v_j.is_boundary())
		{
			int j = v_j.index();
			u_sum += lambda_ij(i,j) * new_positions[j][0];
			v_sum += lambda_ij(i,j) * new_positions[j][1];
		}
	} while(ring.advance());

	Ubar(i) = u_sum;
	Vbar(i) = v_sum;
}

//
// compute mean-value weights given vertices with index i and j
//
float Mesh_fixed_param::lambda_ij(const int i, const int j)
{
 	Mesh_connectivity::Vertex_iterator v_i = mesh().vertex_at(i);
 	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
	float w_sum = 0.0;
	float w_ij;

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
	
		w_ik = (tan(a_ik/2) + tan(b_ki/2)) / r_ik;

		if (k==j) { w_ij = w_ik; }
		w_sum += w_ik;

	} while(ring.advance());

	return w_ij / w_sum;
}

void Mesh_fixed_param::compute_angles(int r_id, Eigen::Vector3d i_pos, Eigen::Vector3d k_pos, float& a_ik, float& b_ki)
{
	Mesh_connectivity::Half_edge_iterator r = mesh().half_edge_at(r_id); // points to i, origin at k

	Eigen::Vector3d p_pos = r.prev().origin().xyz();
	Eigen::Vector3d q_pos = r.next().dest().xyz();

   	float R = (k_pos - i_pos).norm();
	float A1 = (p_pos - i_pos).norm();
	float A2 = (p_pos - k_pos).norm();
	float B1 = (q_pos - i_pos).norm();
	float B2 = (q_pos - k_pos).norm();

	a_ik = acos((R*R + A1*A1 - A2*A2)/(2*R*A1));
	b_ki = acos((B1*B1 + R*R - B2*B2)/(2*R*B1));
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
	
	int n = mesh().n_total_vertices() - boundary.size();
	A = Eigen::MatrixXd::Zero(n, n);
	Ubar = Eigen::MatrixXd::Zero(n, 1);
	Vbar = Eigen::MatrixXd::Zero(n, 1);

    for (const auto& pair : boundary) {
		compute_circle_pos(pair.second);
    }

	// TODO - slow to iterate over all this, maybe a faster solution?
	for(int vid = 0 ; vid < mesh().n_total_vertices() ; ++vid)
	{
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(vid);

		if (!v.is_boundary())
		{
			compute_A_i(vid);
			compute_UVbar_i(vid);
		}
	}
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
