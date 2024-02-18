#include <minimesh/core/mohe/mesh_simplify.hpp>
#include <minimesh/core/util/assert.hpp>
#include <Eigen/Dense> 
#include <set>

namespace minimesh
{
namespace mohe
{

//
// compute Q matrices for all vertices
//
void Mesh_simplify::initialize_Q_matrices()
{
	// printf("q \n");

	for(int v_id = 0 ; v_id < mesh().n_total_vertices(); ++v_id) {
		// printf("q loop\n");
		
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(v_id);
		Mesh_connectivity::Vertex_ring_iterator ring_iter = mesh().vertex_ring_at(v_id);
		Eigen::Matrix4d Q_v = Eigen::Matrix4d::Zero(); 
		
		while(ring_iter.advance()) 
		{
			Mesh_connectivity::Vertex_iterator v_i = ring_iter.half_edge().origin();
			Mesh_connectivity::Vertex_iterator v_j = ring_iter.half_edge().prev().origin();

			Eigen::Vector3d a = v.xyz() - v_i.xyz();
			Eigen::Vector3d b = v_j.xyz() - v_i.xyz();
			Eigen::Vector3d normal = a.cross(b).normalized();
			double A = normal.x();
			double B = normal.y();
			double C = normal.z();
			double D = -normal.dot(v.xyz());

			Eigen::Matrix4d K_v;
			K_v <<
				A * A, A * B, A * C, A * D,
				A * B, B * B, B * C, B * D,
				A * C, B * C, C * C, C * D,
				A * D, B * D, C * D, D * D;

			Q_v += K_v;
		}

		Q_matrices[v_id] = Q_v;
	}
}

//
// compute the new vertex position for an edge to be collapsed, and its associated 
//
void Mesh_simplify::compute_position_and_error(Mesh_connectivity::Half_edge_iterator he)
{
	// printf("pos \n");

	Mesh_connectivity::Vertex_iterator v_1 = he.origin();
	Mesh_connectivity::Vertex_iterator v_2 = he.twin().origin();

	// get Q_v matrix for new vertex 
	Eigen::Matrix4d Q_v = Q_matrices[v_1.index()] + Q_matrices[v_2.index()];
	
	// Assign values to A matrix and b vector from Q_v
	Eigen::Matrix4d A;
	Eigen::Vector4d b;
	A <<
		Q_v(0, 0), Q_v(0, 1), Q_v(0, 2), Q_v(0, 3),
		Q_v(1, 0), Q_v(1, 1), Q_v(1, 2), Q_v(1, 3),
		Q_v(2, 0), Q_v(2, 1), Q_v(2, 2), Q_v(2, 3),
		0.0, 0.0, 0.0, 1.0;
	
	b <<
		0.0, 0.0, 0.0, 1.0;

	// solve system to comput position
	Eigen::Vector4d x = A.fullPivLu().solve(b);

	if (x.allFinite()) 
	{
		new_pos[he.index()] << x(0), x(1), x(2);
	} else {
		// no solution, use midpoint
		new_pos[he.index()] = (v_1.xyz() + v_2.xyz()) / 2;
	}

	double error = x.transpose() * Q_v * x;
	errorQueue.push(Error(error, he.index()));

}


//
// collapse an edge associated with a new vertex 
//
void Mesh_simplify::collapse_edge(Mesh_connectivity::Half_edge_iterator he)
{
	// to be activated
	Mesh_connectivity::Vertex_iterator v = he.origin();
	v.data().half_edge = he.prev().twin().index();
	v.data().xyz = new_pos[he.index()];
	Q_matrices[v.index()] += Q_matrices[he.twin().origin().index()];

	// update origins
	he.next().twin().next().data().origin = v.index();
	he.next().twin().next().twin().next().data().origin = v.index();
	he.twin().prev().twin().prev().twin().data().origin = v.index();
	// Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(he.twin().origin().index());
	// while(ring.advance()) // vertex iterator (gets he pointing TO the vertex)
	// {
	// 	ring.half_edge().twin().data().origin = v.index();
	// }

	// update faces
	he.prev().data().face = he.next().twin().face().index();
	he.prev().face().data().half_edge = he.prev().index();

	he.twin().next().data().face = he.twin().prev().twin().face().index();
	he.twin().next().face().data().half_edge = he.twin().next().index();

	// update edge adjacency	
	he.prev().data().next = he.next().twin().next().index();
	he.prev().data().prev = he.next().twin().prev().index();
	he.next().twin().next().data().prev = he.prev().index();
	he.next().twin().prev().data().next = he.prev().index();
	
	he.twin().next().data().next = he.twin().prev().twin().next().index();
	he.twin().next().data().prev = he.twin().prev().twin().prev().index();
	he.twin().prev().twin().next().data().prev = he.twin().next().index();
	he.twin().prev().twin().prev().data().next = he.twin().next().index();

	// to be deactivate
	// 1 vertex
	he.twin().origin().deactivate();
	
	// 2 faces
	he.face().deactivate();
	he.twin().face().deactivate();

	// adjacent edges
	he.next().twin().origin().data().half_edge = he.prev().index();
	he.next().twin().deactivate();
	he.next().deactivate();
	he.twin().prev().twin().deactivate();
	
	he.twin().prev().origin().data().half_edge = he.twin().next().twin().index();
	he.twin().prev().deactivate();
	
	// collapsed edge
	he.twin().deactivate();
	he.deactivate();

    force_assert( mesh().check_sanity_slowly() );

	// now update position and errors for new adjacent vertices
	Mesh_connectivity::Vertex_ring_iterator v_ring = mesh().vertex_ring_at(v.index());
	while(v_ring.advance()) // vertex iterator (gets he pointing TO the vertex)
	{
		compute_position_and_error(v_ring.half_edge());
	}
}

//
// check to see if the mesh topology is valid. 
// the endpoints of the edge to be collapsed must not have more than two adjacent vertices in common
// true is valid, false if not
//
bool Mesh_simplify::check_topology(Mesh_connectivity::Half_edge_iterator he)
{
	// Create sets to store the indices encountered in each loop
	std::set<int> v1_adj_vs;
	std::set<int> v2_adj_vs;

	// First loop
	Mesh_connectivity::Vertex_ring_iterator v1_ring = mesh().vertex_ring_at(he.origin().index());
	while (v1_ring.advance()) {
		v1_adj_vs.insert(v1_ring.half_edge().origin().index());
	}

	// Second loop
	Mesh_connectivity::Vertex_ring_iterator v2_ring = mesh().vertex_ring_at(he.twin().origin().index());
	while (v2_ring.advance()) {
		v2_adj_vs.insert(v2_ring.half_edge().origin().index());
	}

	// Calculate the intersection of the sets
	std::set<int> intersection;
	std::set_intersection(v1_adj_vs.begin(), v1_adj_vs.end(),
						v2_adj_vs.begin(), v2_adj_vs.end(),
						std::inserter(intersection, intersection.begin()));

	// if interesection != 2, not manifold
	return intersection.size() == 2;
}


// bool Mesh_simplify::check_topology(Mesh_connectivity::Half_edge_iterator he)
// {
// 	// Create sets to store the indices encountered in each loop
// 	std::set<int> v1_adj_vs;
// 	std::set<int> v2_adj_vs;

// 	// First loop
// 	Mesh_connectivity::Vertex_ring_iterator v1_ring = mesh().vertex_ring_at(he.origin().index());
// 	while (v1_ring.advance()) { // vertex iterator (gets he pointing TO the vertex)
// 		Mesh_connectivity::Vertex_iterator vn = v1_ring.half_edge().twin().origin();
// 		Mesh_connectivity::Vertex_ring_iterator vn_ring = mesh().vertex_ring_at(vn.index());
// 		int val;
// 		while (vn_ring.advance()) { // vertex iterator (gets he pointing TO the vertex)
// 			val++;
// 		}
// 		if (val <= 3) { return false; }
// 	}

// 	// Second loop
// 	Mesh_connectivity::Vertex_ring_iterator v2_ring = mesh().vertex_ring_at(he.twin().origin().index());
// 	while (v2_ring.advance()) {
// 		Mesh_connectivity::Vertex_iterator vn = v2_ring.half_edge().twin().origin();
// 		Mesh_connectivity::Vertex_ring_iterator vn_ring = mesh().vertex_ring_at(vn.index());
// 		int val;
// 		while (vn_ring.advance()) { // vertex iterator (gets he pointing TO the vertex)
// 			val++;
// 		}
// 		if (val <= 3) { return false; }
// 	}
// 	return true;
// }

//
// inititialize starting new positions and new errors
//
void Mesh_simplify::init_pos_and_errors() 
{
	printf("init \n");

	for(int he_id = 0 ; he_id < mesh().n_total_half_edges() ; ++he_id) 
	{
		Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);

		if (!he.is_split()) 
		{
			mark_as_split(he);
			compute_position_and_error(he);
		}
	}
	printf("init done \n");

	reset_flags();
}


//
// simplify 
//
void Mesh_simplify::simplify(int num_entities_to_simplify)
{
	printf("qmat size %lu \n", Q_matrices.size());
	for(int entity = 0; entity < num_entities_to_simplify; ++entity) {

		while (!errorQueue.empty()) {

			// printf("not empty \n");
			Error min_error = errorQueue.top();
			Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(min_error.v_id);
			errorQueue.pop();

			// Check if topology is valid
			if (check_topology(he)) {
				// If topology is valid, collapse the edge and exit the loop
				// printf("valid topo \n");

				collapse_edge(he);
				printf("collapsed an edge \n");
				break;
			}
			printf("invalid topo \n");

		}
	}
}

////////////////////////////////////////////////
//////////////////DATA HANDLING/////////////////
////////////////////////////////////////////////

//
// reset any flags
//
void Mesh_simplify::reset_flags() 
{
	// reset split edge flag
	while(!split_half_edges.empty()) {
		int index = split_half_edges.top();
		Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(index);
		he.data().is_split = false;
		split_half_edges.pop();
	}
	assert(split_half_edges.empty());

	// TODO MAYBE NEED TO DELETE ERRORQUEUE?

}

//
// marks given half edge and its twin as split
//
void Mesh_simplify::mark_as_split(Mesh_connectivity::Half_edge_iterator he)
{
	he.data().is_split = true;
	split_half_edges.push(he.index());

	he.twin().data().is_split = true;
	split_half_edges.push(he.twin().index());
}

} // end of mohe
} // end of minimesh
