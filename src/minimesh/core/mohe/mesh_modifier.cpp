#include <minimesh/core/mohe/mesh_modifier.hpp>
#include <minimesh/core/util/assert.hpp>

namespace minimesh
{
namespace mohe
{

//
// Given a half edge, split the edge into two
//
void Mesh_modifier::edge_split(const int he_index)
{
	// get he from he_index
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_index);
	Mesh_connectivity::Half_edge_iterator he_twin = he.twin();
	Mesh_connectivity::Half_edge_iterator he_twin_next = he.twin().next();
	Mesh_connectivity::Half_edge_iterator he_prev = he.prev();
	Mesh_connectivity::Half_edge_iterator he_next = he.next();

	// add a new vertex v to mesh
	Mesh_connectivity::Vertex_iterator v = mesh().add_vertex();
	v.data().half_edge = he.index();
	// v.data().xyz = (he.origin().data().xyz + he_twin.origin().data().xyz) / 2; 
	v.data().is_new = true;
	new_vertices.push(v.index());

	// create two new half edges for split
	Mesh_connectivity::Half_edge_iterator he_new = mesh().add_half_edge();
	Mesh_connectivity::Half_edge_iterator he_new_twin = mesh().add_half_edge();

	Mesh_connectivity::Vertex_iterator v_old = mesh().vertex_at(he.data().origin);
	v_old.data().half_edge = he_new.index();

	// update he_new
	he_new.data().face = he.data().face;
	he_new.data().next = he.index();
	he_new.data().prev = he_prev.index();
	he_new.data().twin = he_new_twin.index();
	he_new.data().origin = he.data().origin;
	he_new.data().is_split = true;
	split_half_edges.push(he_new.index());

	// update he: face, next, twin unchanged
	he.data().prev = he_new.index();
	he.data().origin = v.index();
	he.data().is_split = true;
	split_half_edges.push(he.index());

	// update he_prev: all but next unchanged
	he_prev.data().next = he_new.index();

	// update he_new_twin
	he_new_twin.data().face = he_twin.data().face;
	he_new_twin.data().next = he_twin.data().next;
	he_new_twin.data().prev = he_twin.index();
	he_new_twin.data().twin = he_new.index();
	he_new_twin.data().origin = v.index();
	he_new_twin.data().is_split = true;
	split_half_edges.push(he_new_twin.index());

	// udpate he_twin: face, prev, twin, origin unchanged
	he_twin.data().next = he_new_twin.index();
	he_twin.data().is_split = true;
	split_half_edges.push(he_twin.index());

	// update he_twin_next: all but prev unchanged
	he_twin_next.data().prev = he_new_twin.index();
}

//
// create additional faces. cut the first corner found, unless the face is already a triangle
//
void Mesh_modifier::cut_a_corner(const int face_index)
{
	// get face from face_index
	Mesh_connectivity::Face_iterator face = mesh().face_at(face_index);

	// get half edge of face
	Mesh_connectivity::Half_edge_iterator he = face.half_edge();

	// make sure half edge points at a corner
	while(he.twin().origin().is_new()) {
 		he = he.next();
		// printf("next");
	}

	// get half edge neighbour data
	Mesh_connectivity::Half_edge_iterator he_prev = mesh().half_edge_at(he.data().prev);
	Mesh_connectivity::Half_edge_iterator he_next = mesh().half_edge_at(he.data().next);
	Mesh_connectivity::Half_edge_iterator he_next_next = mesh().half_edge_at(he_next.data().next);

	// add he_new to mesh 
	Mesh_connectivity::Half_edge_iterator he_new = mesh().add_half_edge();
	Mesh_connectivity::Half_edge_iterator he_new_twin = mesh().add_half_edge();

	// create new face for corner
	Mesh_connectivity::Face_iterator face_new = mesh().add_face();
	face_new.data().half_edge = he_new.index();

	// update he_new
	he_new.data().face = face_new.index();
	he_new.data().next = he.index();
	he_new.data().prev = he_next.index();
	he_new.data().twin = he_new_twin.index();
	he_new.data().origin = he_next_next.data().origin;
	he_new.data().is_split = true;
	split_half_edges.push(he_new.index());

	// update he
	he.data().face = face_new.index();
	he.data().prev = he_new.index();

	// update he_next
	he_next.data().face = face_new.index();
	he_next.data().next = he_new.index();

	// update he_new_twin
	he_new_twin.data().face = face_index;
	he_new_twin.data().next = he_next_next.index();
	he_new_twin.data().prev = he_prev.index();
	he_new_twin.data().twin = he_new.index();
	he_new_twin.data().origin = he.data().origin;
	he_new_twin.data().is_split = true;
	split_half_edges.push(he_new_twin.index());

	// update he_prev
	he_prev.data().next = he_new_twin.index();

	// update he_next_next 
	he_next_next.data().prev = he_new_twin.index();

	// update orig face data
	face.data().half_edge = he_new_twin.index();
}

//
// checks to see if a given face is a triangle
//
bool Mesh_modifier::is_triangle(const int face_index)
{
	// get face from face_index
	Mesh_connectivity::Face_iterator face = mesh().face_at(face_index);

	// get index of half edge of face
	int he_index = face.data().half_edge;
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_index);
	Mesh_connectivity::Half_edge_iterator he_next = mesh().half_edge_at(he.data().next);
	Mesh_connectivity::Half_edge_iterator he_next_next = mesh().half_edge_at(he_next.data().next);

	// check if face is minimal triangle
	return he.data().prev == he_next.data().next;
}

//
// butterfly schema mask for new vertex pos
//
Eigen::Vector3d Mesh_modifier::butterfly_mask(const int he_index)
{
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_index);
	he.data().is_split = true;
	split_half_edges.push(he.index());

	he.twin().data().is_split = true;
	split_half_edges.push(he.twin().index());

	Eigen::Vector3d ml = he.origin().xyz(); 
	Eigen::Vector3d mr = he.twin().origin().xyz(); 

	Eigen::Vector3d tl = he.prev().twin().prev().origin().xyz(); 
	Eigen::Vector3d tm = he.prev().origin().xyz(); 
	Eigen::Vector3d tr = he.next().twin().prev().origin().xyz(); 
	
	Eigen::Vector3d bl = he.twin().next().twin().prev().origin().xyz();
	Eigen::Vector3d bm = he.twin().prev().origin().xyz(); 
	Eigen::Vector3d br = he.twin().prev().twin().prev().origin().xyz(); 
	
	return (-1.0/16.0) * tl + (2.0/16.0) * tm + (-1.0/16.0) * tr +
		(8.0/16.0) * ml + (8.0/16.0) * mr +
		(-1.0/16.0) * bl + (2.0/16.0) * bm + (-1.0/16.0) * br;
} 

//
// reset the flags applied in subdivision step. eg. is_new and is_split
//
void Mesh_modifier::reset_flags() 
{
	// reset split edge flag
	while(!split_half_edges.empty()) {
		int index = split_half_edges.top();
		Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(index);
		he.data().is_split = false;
		split_half_edges.pop();
	}

}

//
// subdivision 
//
void Mesh_modifier::subdivision()
{

	// compute new vertex pos using butterfly stencil mask
	std::queue<Eigen::Vector3d> new_vertices_pos;
	for(int he_id = 0 ; he_id < mesh().n_total_half_edges() ; ++he_id)
	{
		Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);

		// prevent computing new vertex pos for twins
		if (he.is_active() && !he.is_split())
		{
			new_vertices_pos.push(butterfly_mask(he_id));
		}
	}

	reset_flags();

	printf("splitting edges \n");
	// iterate over all half edges that are not split
	for(int he_id = 0 ; he_id < mesh().n_total_half_edges() ; ++he_id)
	{
		Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);

		// split half edge if not split
		if (he.is_active() && !he.is_split())
		{
			edge_split(he_id);
		}
	}

    force_assert( mesh().check_sanity_slowly() );

	printf("cutting corners \n");
	// iterate over all non-triangular faces
	for(int face_id = 0 ; face_id < mesh().n_total_faces() ; ++face_id)
	{
		Mesh_connectivity::Face_iterator face = mesh().face_at(face_id);

		// cut a corner out of the face if it is not a traingle
		while(face.is_active() && !is_triangle(face_id))
		{
			cut_a_corner(face_id);
		}
	}

    force_assert( mesh().check_sanity_slowly() );

	// update new vertices pos
	// index = 0;
	// for(int v_id = 0 ; v_id < mesh().n_total_vertices() ; ++v_id)
	// {
	// 	Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(v_id);
	// 	if(v.is_active() && v.is_new())
	// 	{
	// 		v.data().xyz = new_vertex_pos[index];
	// 		++index;
	// 	}
	// }

	// printf("new_vertex_queue_size = %d \n", new_vertices.size());
	// printf("new_pos_queue_size = %d \n", new_vertices_pos.size());


	while(!new_vertices.empty()) {
		int vid = new_vertices.front();
		Eigen::Vector3d pos = new_vertices_pos.front();

		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(vid);
		v.data().is_new = false;
		v.data().xyz = pos;

		new_vertices.pop();
		new_vertices_pos.pop();
	}



	reset_flags();
}

//
// Given two vertices, this function return the index of the half-edge going from v0 to v1.
// Returns -1 if no half-edge exists between the two vertices.
//
int Mesh_modifier::get_halfedge_between_vertices(const int v0, const int v1)
{
	// Get a ring iterator for v0
	Mesh_connectivity::Vertex_ring_iterator ring_iter = mesh().vertex_ring_at(v0);

	int answer = mesh().invalid_index;

	// Loop over all half-edges that end at v0.
	do
	{
		// Make sure that the half-edge does end at v0
		assert(ring_iter.half_edge().dest().index() == v0);

		// If the half-edge also starts and v1, then it's twin
		// goes from v0 to v1. This would be the half-edge that
		// we were looking for
		if(ring_iter.half_edge().origin().index() == v1)
		{
			answer = ring_iter.half_edge().twin().index();
		}
	} while(ring_iter.advance());

	if(answer != mesh().invalid_index)
	{
		assert(mesh().half_edge_at(answer).origin().index() == v0);
		assert(mesh().half_edge_at(answer).dest().index() == v1);
	}

	return answer;
}


bool Mesh_modifier::flip_edge(const int he_index)
{
	//
	// Take a reference to all involved entities
	//

	// HALF-EDGES
	Mesh_connectivity::Half_edge_iterator he0 = mesh().half_edge_at(he_index);
	Mesh_connectivity::Half_edge_iterator he1 = he0.twin();

	// meshes on the boundary are not flippable
	if(he0.face().is_equal(mesh().hole()) || he1.face().is_equal(mesh().hole()))
	{
		return false;
	}

	Mesh_connectivity::Half_edge_iterator he2 = he0.next();
	Mesh_connectivity::Half_edge_iterator he3 = he2.next();
	Mesh_connectivity::Half_edge_iterator he4 = he1.next();
	Mesh_connectivity::Half_edge_iterator he5 = he4.next();

	// VERTICES
	Mesh_connectivity::Vertex_iterator v0 = he1.origin();
	Mesh_connectivity::Vertex_iterator v1 = he0.origin();
	Mesh_connectivity::Vertex_iterator v2 = he3.origin();
	Mesh_connectivity::Vertex_iterator v3 = he5.origin();

	// FACES
	Mesh_connectivity::Face_iterator f0 = he0.face();
	Mesh_connectivity::Face_iterator f1 = he1.face();

	//
	// Now modify the connectivity
	//

	// HALF-EDGES
	he0.data().next = he3.index();
	he0.data().prev = he4.index();
	he0.data().origin = v3.index();
	//
	he1.data().next = he5.index();
	he1.data().prev = he2.index();
	he1.data().origin = v2.index();
	//
	he2.data().next = he1.index();
	he2.data().prev = he5.index();
	he2.data().face = f1.index();
	//
	he3.data().next = he4.index();
	he3.data().prev = he0.index();
	//
	he4.data().next = he0.index();
	he4.data().prev = he3.index();
	he4.data().face = f0.index();
	//
	he5.data().next = he2.index();
	he5.data().prev = he1.index();

	// VERTICES
	v0.data().half_edge = he2.index();
	v1.data().half_edge = he4.index();
	v2.data().half_edge = he1.index();
	v3.data().half_edge = he0.index();

	// FACES
	f0.data().half_edge = he0.index();
	f1.data().half_edge = he1.index();

	// operation successful
	return true;
} // All done


} // end of mohe
} // end of minimesh
