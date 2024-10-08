#include <minimesh/core/mohe/mesh_subdivision.hpp>
#include <minimesh/core/util/assert.hpp>
#include "minimesh/core/util/numbers.hpp"

namespace minimesh
{
namespace mohe
{

//
// Given a half edge, split the edge into two
//
void Mesh_subdivision::edge_split(const int he_index)
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
void Mesh_subdivision::cut_a_corner(const int face_index)
{
	// get face from face_index
	Mesh_connectivity::Face_iterator face = mesh().face_at(face_index);

	// get half edge of face
	Mesh_connectivity::Half_edge_iterator he = face.half_edge();

	// make sure half edge points at a corner
	while(he.twin().origin().is_new()) {
 		he = he.next();
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
bool Mesh_subdivision::is_triangle(const int face_index)
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
// checks to see if a given half edge is at a boundary
//
bool Mesh_subdivision::is_boundary(const int he_index)
{
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_index);

	// see if the face for given half edge (or its twin) has negative index
	return he.twin().face().index() < 0 || he.face().index() < 0; 
}

//
// loop schema mask for control interior vertex
//
Eigen::Vector3d Mesh_subdivision::loop_interior_control(const int v_index)
{
	Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(v_index);
	Mesh_connectivity::Vertex_ring_iterator ring_iter = mesh().vertex_ring_at(v_index);
	v.data().is_done = true;

	int n = 1; // valence of vertex
	Eigen::Vector3d sum_adj_pos = ring_iter.half_edge().origin().xyz();
	
	while(ring_iter.advance()) 
	{
		sum_adj_pos += ring_iter.half_edge().origin().xyz();
		n++;
	}

    double alpha_n = (3.0/8.0) + std::pow((3.0/8.0) + (1.0/4.0) * std::cos((2.0 * numbers::pi) / n), 2);

	return alpha_n * v.xyz() + ((1.0-alpha_n)/n) * sum_adj_pos;
} 


//
// loop schema mask for new interior vertex
//
Eigen::Vector3d Mesh_subdivision::loop_interior_new(const int he_index)
{
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_index);
	he.data().is_split = true;
	split_half_edges.push(he.index());

	he.twin().data().is_split = true;
	split_half_edges.push(he.twin().index());

	Eigen::Vector3d ml = he.origin().xyz(); 
	Eigen::Vector3d mr = he.twin().origin().xyz(); 

	Eigen::Vector3d tm = he.prev().origin().xyz(); 
	Eigen::Vector3d bm = he.twin().prev().origin().xyz(); 
	
	return (3.0/8.0) * (ml + mr) + (1.0/8.0) * (tm + bm);
} 

//
// loop schema mask for control boundary vertex
//
Eigen::Vector3d Mesh_subdivision::loop_boundary_control(const int he_index)
{
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_index);

	Mesh_connectivity::Vertex_iterator mid = he.origin();
	mid.data().is_done = true; // so we don't recompute when revisited
	
	// boundary edges are a trivial cycle so we can just use next/prev
	Mesh_connectivity::Half_edge_iterator right = he.next();
	Mesh_connectivity::Half_edge_iterator left = he.prev();

	return (1.0/8.0) * left.origin().xyz() + (6.0/8.0) * mid.xyz() + (1.0/8.0) * right.origin().xyz();
} 


//
// loop schema mask for new boundary vertex 
//
Eigen::Vector3d Mesh_subdivision::loop_boundary_new(const int he_index)
{
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_index);
	he.data().is_split = true;
	split_half_edges.push(he.index());

	he.twin().data().is_split = true;
	split_half_edges.push(he.twin().index());

	Eigen::Vector3d v_left = he.origin().xyz(); 
	Eigen::Vector3d v_right = he.twin().origin().xyz(); 
	
	return (1.0/2.0) * (v_left + v_right);
} 

//
// reset the flags applied in subdivision step. eg. is_split
//
void Mesh_subdivision::reset_flags() 
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

//
// subdivision 
//
void Mesh_subdivision::subdivision()
{
	int n_contol_vertices = mesh().n_total_vertices();

	// precompute new half edge vertex pos with desired stencil mask
	for(int he_id = 0 ; he_id < mesh().n_total_half_edges() ; ++he_id)
	{
		Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);

		// half edge is on boundary
		if (is_boundary(he_id)) 
		{
			// prevent computing new vertex pos for twins
			if (!he.is_split())
			{
				// compute position for new boundary vertex 
				new_vertices_pos.push(loop_boundary_new(he_id));
				
				// get the actual boundary half edge (not the interior twin)
				if (!(he.face().index() < 0)) { he = he.twin();}

				// compute postion for control boundary vertex
				vertex_pos_map[he.origin().index()] = loop_boundary_control(he.index());
			}
		} 
		else 
		{
			// prevent computing new vertex pos for twins
			if (!he.is_split())
			{
				// compute position for new interior vertex 
				new_vertices_pos.push(loop_interior_new(he_id));
				
				// compute position for control vertices (adjacent to edge)
				Mesh_connectivity::Vertex_iterator v1 = he.origin();
				Mesh_connectivity::Vertex_iterator v2 = he.twin().origin();

				// dont bother computing boundary vertex 
				if (!v1.is_done())
				{
					vertex_pos_map[v1.index()] = loop_interior_control(v1.index());
				}

				if (!v2.is_done())
				{
					vertex_pos_map[v2.index()] = loop_interior_control(v2.index());
				}
			}
		}	
	}

	reset_flags();

	// split half edges and add new vertex
	for(int he_id = 0 ; he_id < mesh().n_total_half_edges() ; ++he_id)
	{
		Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);
		if (he.is_active() && !he.is_split())
		{
			edge_split(he_id);
		}
	}

    force_assert( mesh().check_sanity_slowly() );

	// add half edges to cut corners of non-triangular faces
	for(int face_id = 0 ; face_id < mesh().n_total_faces() ; ++face_id)
	{
		Mesh_connectivity::Face_iterator face = mesh().face_at(face_id);
		while(face.is_active() && !is_triangle(face_id))
		{
			cut_a_corner(face_id);
		}
	}

	force_assert( mesh().check_sanity_slowly() );

	// update control vertex positions (also clear flags/data structures)
	for(int v_id = 0 ; v_id < n_contol_vertices ; ++v_id)
	{
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(v_id);
		v.data().is_done = false;
		mesh().vertex_at(v_id).data().xyz = vertex_pos_map[v_id];
		vertex_pos_map.erase(v_id);
	}

	assert(vertex_pos_map.empty());

	// postion new vertices (also clear flags/data structures)
	while(!new_vertices.empty()) {
		int vid = new_vertices.front();
		Eigen::Vector3d pos = new_vertices_pos.front();

		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(vid);
		v.data().is_new = false;
		v.data().xyz = pos;

		new_vertices.pop();
		new_vertices_pos.pop();
	}

	assert(new_vertices.empty());
	assert(new_vertices_pos.empty());

	reset_flags();
}


} // end of mohe
} // end of minimesh
