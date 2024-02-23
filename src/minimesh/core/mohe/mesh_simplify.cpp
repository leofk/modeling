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
	for(int v_id = 0 ; v_id < mesh().n_total_vertices(); ++v_id) {
		
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(v_id);
		Mesh_connectivity::Vertex_ring_iterator ring_iter = mesh().vertex_ring_at(v_id);
		Eigen::Matrix4d Q_v = Eigen::Matrix4d::Zero(); 
		
		do {
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
		} while(ring_iter.advance()); 

		Q_matrices[v_id] = Q_v;
	}
}

//
// compute the new vertex position for an edge to be collapsed, and its associated 
//
void Mesh_simplify::compute_position_and_error(Mesh_connectivity::Half_edge_iterator he)
{
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
		new_pos[he.twin().index()] << x(0), x(1), x(2);
	} else {
		// no solution, use midpoint
		new_pos[he.index()] = (v_1.xyz() + v_2.xyz()) / 2;
		new_pos[he.twin().index()] = (v_1.xyz() + v_2.xyz()) / 2;
	}

	double error = x.transpose() * Q_v * x;
	errorQueue.push(Error(error, he.index()));
	errorMap[he.index()] = error;
}

//
// collapse an edge associated with a new vertex 
//
void Mesh_simplify::collapse_edge(Mesh_connectivity::Half_edge_iterator he)
{
	Mesh_connectivity::Half_edge_iterator H = he;
	Mesh_connectivity::Half_edge_iterator HT = he.twin();
	Mesh_connectivity::Half_edge_iterator HP = he.prev();
	Mesh_connectivity::Half_edge_iterator HN = H.next();
	Mesh_connectivity::Half_edge_iterator HTP = HT.prev();
	Mesh_connectivity::Half_edge_iterator HTN = HT.next();
	Mesh_connectivity::Half_edge_iterator HNT = HN.twin();
	Mesh_connectivity::Half_edge_iterator HTPT = HTP.twin();
	Mesh_connectivity::Half_edge_iterator HNTP = HNT.prev();
	Mesh_connectivity::Half_edge_iterator HNTN = HNT.next();
	Mesh_connectivity::Half_edge_iterator HTPTN = HTPT.next();
	Mesh_connectivity::Half_edge_iterator HTPTP = HTPT.prev();

	Mesh_connectivity::Vertex_iterator V1 = H.origin();
	Mesh_connectivity::Vertex_iterator V2 = H.dest();
	Mesh_connectivity::Vertex_iterator V3 = HN.dest();
	Mesh_connectivity::Vertex_iterator V4 = HTP.origin();

	Mesh_connectivity::Face_iterator F1 = H.face();
	Mesh_connectivity::Face_iterator F2 = HT.face();
	Mesh_connectivity::Face_iterator F3 = HNT.face();
	Mesh_connectivity::Face_iterator F4 = HTPT.face();

	// check H was not deactivated in a previous iteration
	if (!H.is_active()) return;

	// check no vertices adjacent to H has 3 or less valence
	if (!check_valence(H)) return;

	// check that collapsing an edge wont result in illegal topology
	if (!check_connectivity(H)) return;

	// topology will remain valid, now continue with mesh modifications

	// update origins
	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(V2.index());
	do // *points to*
	{
		ring.half_edge().twin().data().origin = V1.index();
	} while(ring.advance());

	// update edge adjacency (next/prevs)
	HP.data().next = HNTN.index();
	HP.data().prev = HNTP.index();
	HNTN.data().prev = HP.index();
	HNTP.data().next = HP.index();

	HTN.data().next = HTPTN.index();
	HTN.data().prev = HTPTP.index();
	HTPTN.data().prev = HTN.index();
	HTPTP.data().next = HTN.index();

	// update vertex edge adjacency
	V3.data().half_edge = HP.index(); // *going out*
	V4.data().half_edge = HTPTN.index(); // *going out*

	// update faces
	HP.data().face = F3.index();
	F3.data().half_edge = HP.index();

	HTN.data().face = F4.index();
	F4.data().half_edge = HTN.index();

	// update V1 
	V1.data().half_edge = HNTN.index(); // *going out*
	V1.data().xyz = new_pos[he.index()];
	Q_matrices[V1.index()] += Q_matrices[V2.index()];

	// deactivate 1 vertex
	V2.deactivate();
	
	// deactivate 2 faces
	F1.deactivate();
	F2.deactivate();

	// deactivate 6 half-edges
	H.deactivate();
	HT.deactivate();
	HN.deactivate();
	HTP.deactivate();
	HNT.deactivate();
	HTPT.deactivate();

    // force_assert( mesh().check_sanity_slowly() );

	// now update position and errors for new adjacent vertices
	ring = mesh().vertex_ring_at(V1.index());
	do // vertex iterator (gets he pointing TO the vertex)
	{
		compute_position_and_error(ring.half_edge());
	} while(ring.advance());
}

//
// only one edge should merge after a collapse
// see that end point vertices only have two adjacent vertices in common
// true if connectivity is valid
//
bool Mesh_simplify::check_connectivity(Mesh_connectivity::Half_edge_iterator he)
{
	// Create sets to store the indices encountered in each loop
	std::set<int> v1_adj_vs;
	std::set<int> v2_adj_vs;

	// First loop
	Mesh_connectivity::Vertex_ring_iterator v1_ring = mesh().vertex_ring_at(he.origin().index());
	do // *points to*
	{
		v1_adj_vs.insert(v1_ring.half_edge().origin().index());
	} while(v1_ring.advance()); 

	// Second loop
	Mesh_connectivity::Vertex_ring_iterator v2_ring = mesh().vertex_ring_at(he.twin().origin().index());
	do
	{
		v2_adj_vs.insert(v2_ring.half_edge().origin().index());
	} while(v2_ring.advance()); 


	// Calculate the intersection of the sets
	std::set<int> intersection;
	std::set_intersection(v1_adj_vs.begin(), v1_adj_vs.end(),
						v2_adj_vs.begin(), v2_adj_vs.end(),
						std::inserter(intersection, intersection.begin()));

	// if interesection != 2, not manifold
	return intersection.size() == 2;
}

//
// valence of endpoint vertices of edge being collapsed should not be 3 
// true if valence > 3
//
bool Mesh_simplify::check_valence(Mesh_connectivity::Half_edge_iterator he)
{
	int val_v2 = 0;
	Mesh_connectivity::Vertex_ring_iterator v2_ring = mesh().vertex_ring_at(he.dest().index());
	do
	{
		val_v2++;
	} while(v2_ring.advance()); 

	return val_v2 > 3;
}


//
// inititialize starting new positions and new errors
//
void Mesh_simplify::init_pos_and_errors() 
{
	for(int he_id = 0 ; he_id < mesh().n_total_half_edges() ; ++he_id) 
	{
		Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);

		if (!he.is_split()) 
		{
			mark_as_split(he);
			compute_position_and_error(he);
		}
	}

}

Eigen::Vector4f generateColorFromSpectrum(float position, float maxPosition) {
    float hue = (position / maxPosition) * 240.0f; // Map position to hue value (0 to 240 degrees)
    float saturation = 1.0f; // Full saturation
    float value = 1.0f; // Full value (brightness)

    // Convert HSV to RGB
    int hi = static_cast<int>(std::floor(hue / 60.0f)) % 6;
    float f = hue / 60.0f - std::floor(hue / 60.0f);
    float v = value;
    float p = value * (1 - saturation);
    float q = value * (1 - f * saturation);
    float t = value * (1 - (1 - f) * saturation);

    switch (hi) {
        case 0:
            return Eigen::Vector4f(v, t, p, 1.0f);
        case 1:
            return Eigen::Vector4f(q, v, p, 1.0f);
        case 2:
            return Eigen::Vector4f(p, v, t, 1.0f);
        case 3:
            return Eigen::Vector4f(p, q, v, 1.0f);
        case 4:
            return Eigen::Vector4f(t, p, v, 1.0f);
        default:
            return Eigen::Vector4f(v, p, q, 1.0f);
    }
}


void Mesh_simplify::color_queue(Mesh_buffer &buffer, Mesh_connectivity::Defragmentation_maps &defrag, int num_e) 
{
	if (errorQueue.empty()) return;
    std::priority_queue<Error> temp = errorQueue;
	Eigen::Matrix4Xf vertexColors = Eigen::Matrix4Xf::Ones(4, mesh().n_active_vertices());
	int max = std::min((int) errorQueue.size(), std::min(30, num_e));
	for (int i = max; i > 0; --i) {
		float position = static_cast<float>(max - i + 1); 
        Eigen::Vector4f rgba = generateColorFromSpectrum(position, max); 
		Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(temp.top().he_id);

		int v1 = defrag.old2new_vertices[he.origin().index()];
		int v2 = defrag.old2new_vertices[he.dest().index()];
		if (v1 != mesh().invalid_index) vertexColors.col(v1) << rgba[0], rgba[1], rgba[2], rgba[3];
		if (v2 != mesh().invalid_index) vertexColors.col(v2) << rgba[0], rgba[1], rgba[2], rgba[3];

		temp.pop();
	}

	buffer.set_vertex_colors(vertexColors);
}


//
// simplify 
//
void Mesh_simplify::simplify(int num_entities_to_simplify)
{
	for(int entity = 0; entity < num_entities_to_simplify; ++entity) {
		while (!errorQueue.empty()) {

			Error min_error = errorQueue.top();
			Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(min_error.he_id);
			errorQueue.pop();
			
			// make sure the error is the most recent for the given edge
			if (min_error.value == errorMap[he.index()]) {
				collapse_edge(he);
				break;
			}
		}
		if(errorQueue.empty()) 
		{
			reset_flags();
			break;
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

    // Clear Q_matrices
    Q_matrices.clear();

    // Clear new_pos
    new_pos.clear();

    // Clear errorMap
    errorMap.clear();
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
