#include <minimesh/core/mohe/mesh_mapping.hpp>
#include <minimesh/core/util/assert.hpp>

namespace minimesh
{
namespace mohe
{

void Mesh_mapping::build_mapping()
{
	// inter surface mapping

	// given two meshes, one-to-one vertex correspondance, same connectivity'
	// m1 - mesh 1 as it were
	// m1p - mesh 2 on mesh 1s positions
	// m2 - mesh 2 as it were
	// m2p - mesh 1 on mesh 2s positions

	// init
		// maps from m1 -> m2p
		// maps from m2 -> m1p
		// intersurface map

	// todo - init
	/*
	- we split a vertex on m1, call it v1
	- we use the same info to split vertex on m2'
	- now we need to figure out where m2' is on m2 to derive the inter-surface mapping
	- now split a vertex on m2
	- we need to check if this split affected any mapping from m1->m2, i.e if faces changed
	- go to m1' and split vert
	- compute mapping for vert on m1' is relative to m1
	- repeat
	*/
}

//
// split the vertex in original mesh according to reverse edge-collapse
//
// void Mesh_mapping::vertex_split()
// {
	// suppose we split on m1
	// here we need to check if this split fucked with our previous mappings from m2-> m1
		// ie. suppose v on m2, there is a map v -> some face in m1, but b/c of the 
		// split on m1, v might no longer belong to that face. so we need to realize its bary cords
	// now compute barycentric coordinates of new v (in m1) with respect to its ring
	// use map to find positions N(v) 
	// compute v' position on m2 using bary ^^ maybe above can be parrallel	
	// v' must be contained in of the faces (in m2) in map used to find N(v)
	// see which of the faces contain v' (via bary coords)
	// the face with non negative coords is our mapping
// }

void Mesh_mapping::ISM_iteration(int m)
{
	// say M1, split vertex. 
	std::map<int, inter_data> & ISM = ISM_M1;
	if (m == 1) {
		Mesh_connectivity & mesh = _m1;
	}
	if (m == 2) {
		Mesh_connectivity & mesh = _m2;
		ISM = ISM_M2;
	}

	// 1. pop from top of collapse queue - HistroyEntry

	//todo adapt akumah

	// 2. get keptFaces and removedVertex from history Entry
	std::vector<int> faces; // = historyEntry.keptFaces
	Eigen::Vector3d v_pos;  // = historyEntry.removedVertex.position
	int v_id;  // = historyEntry.removedVertex.id

	// 3. compute bary coords on keptFaces vs. removedVertex to see where removedVertex lives
	inter_data v_data = find_enclosing_face(faces, v_pos);

	// 4. once we find the face, get its vertices, and use the map to find their positions in m2'
	int f_id = v_data.f_id;
	Mesh_connectivity::Face_iterator f = mesh().face_at(f_id);
	Mesh_connectivity::Half_edge_iterator he = f.half_edge();
	Mesh_connectivity::Vertex_iterator v = he.origin();
	std::vector<int> Mj_faces;
	std::vector<Eigen::Vector3d> Mjp_verts;

	do
	{
		Mesh_connectivity::Vertex_iterator v_curr = he.origin();
		inter_data curr_data = ISM[v_curr.index()];
		Mj_faces.push_back(curr_data.f_id);
		Mjp_verts.push_back(get_pos_from_inter_data(curr_data));
		he = he.next();
	} while (he.origin().index() != v.index());

	// 5. compute v' using bary coords C and positions in m2'
	Eigen::Vector3d mj_vp;
	for (int i=0; i<3; i++)
	{
		mj_vp += Mjp_verts[i]*v_data.bary[i];
	}

	// now we need to figure out what face in m2 v' is in
	// 6. collect the faces on m2 from the map used in 4
	// 7. for each distinct face from 6. lets compute bary coords for v' 
	//	  if the coords are non-negative, this is our face F.
	inter_data vp_data = find_enclosing_face(Mj_faces, mj_vp);

	// 8. add v to the map m1->m2' or v1 -> F, C
	ISM[v_id] = vp_data;
}

//
// compute position of vertice using barycentric coordinates given by inter data
//
Eigen::Vector3d Mesh_mapping::get_pos_from_inter_data(inter_data data)
{
	int f_id = data.f_id;
	Mesh_connectivity::Face_iterator f = mesh().face_at(f_id);
	Mesh_connectivity::Half_edge_iterator he = f.half_edge();
	Mesh_connectivity::Vertex_iterator v = he.origin();
	std::vector<Eigen::Vector3d> verts;

	do
	{
		Mesh_connectivity::Vertex_iterator v_curr = he.origin();
		verts.push_back(v_curr.xyz());
		he = he.next();
	} while (he.origin().index() != v.index());

	Eigen::Vector3d v_out;
	for (int i=0; i<3; i++)
	{
		v_out += verts[i]*data.bary[i];
	}
	return v_out;
}

//
// TODO - compute the barycentric coordinates of v_pos wihin each face in faces
// return the face and coordinates of the face with non-negative bary. coords
//
inter_data Mesh_mapping::find_enclosing_face(std::vector<int> faces, Eigen::Vector3d v_pos)
{
	std::vector<double> bary; 
	int f_id;
	return inter_data(f_id, bary);
}


void Mesh_mapping::init_ISM()
{
	for(auto &pair: mesh_map){
		ISM_M1[pair.first] = trivial_map_data(pair.second);
		ISM_M2[pair.second] = trivial_map_data(pair.first);
	}
}

inter_data Mesh_mapping::trivial_map_data(int vid)
{
	std::vector<double> bary; 

	Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(vid);
	Mesh_connectivity::Half_edge_iterator he_v = v.half_edge();
	Mesh_connectivity::Face_iterator f_v = he_v.face();

	Mesh_connectivity::Half_edge_iterator f_he = f_v.half_edge();
	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(f_he.origin().index());

	do
	{
		if (ring.half_edge().origin().index() == vid) bary.push_back(1.0);
		else bary.push_back(0.0);
	} while (ring.advance());
	
	int f_id = f_v.index(); 
	return inter_data(f_id, bary);
}

//
// barycentric coordinates (mean-value)
// the weight of vertex j for i
//
// double Mesh_fixed_param::get_wik(int he_index) 
// {
// 	// points to i from j
// 	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_index); 

// 	Mesh_connectivity::Vertex_iterator I = he.dest();
// 	Mesh_connectivity::Vertex_iterator J = he.origin();  
// 	Mesh_connectivity::Vertex_iterator P = he.twin().next().dest(); 
// 	Mesh_connectivity::Vertex_iterator Q = he.next().dest(); 

// 	Eigen::Vector3d I_xyz = I.xyz();
// 	Eigen::Vector3d J_xyz = J.xyz();
// 	Eigen::Vector3d P_xyz = P.xyz();
// 	Eigen::Vector3d Q_xyz = Q.xyz();

// 	Eigen::Vector3d IJ = J_xyz-I_xyz;
// 	Eigen::Vector3d IP = P_xyz-I_xyz;
// 	Eigen::Vector3d IQ = Q_xyz-I_xyz;

// 	double alpha = get_angle(IJ, IP);
// 	double beta = get_angle(IJ, IQ);
// 	double r = IJ.norm();

// 	return (std::tan(alpha / 2) + std::tan(beta / 2)) / r;
// }

} // end of mohe
} // end of minimesh
