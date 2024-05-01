#include <minimesh/core/mohe/mesh_mapping.hpp>
#include <minimesh/core/mohe/mesh_simplify.hpp>
#include <minimesh/core/util/assert.hpp>
#include <set>
#include <iostream>

namespace minimesh
{
namespace mohe
{

//
// Build Inter-surface map between two meshes
//
void Mesh_mapping::build_mapping()
{
	// Simplify Cow to 4 verts
    Mesh_simplify simp_m1(_m1);
	simp_m1.init();
	simp_m1.simplify_to_target(4);
	printf("M1 Simplified. \n");

	// Simplify Camel to 4 verts
    Mesh_simplify simp_m2(_m2);
	simp_m2.init();
	simp_m2.simplify_to_target(4);
	printf("M2 Simplified. \n");

	// Initialize Inter-surface Map with initial trivial mapping
	init_ISM();

	printf("Begin Iteration. \n");
	int count = 0;
	int breakpoint = 3;

	// Iteratively reconstruct both meshes 
	while (!simp_m1.is_history_empty() || !simp_m2.is_history_empty()) 
	{
		if (!simp_m1.is_history_empty()) {
			current_mesh = M1;
			split_vertex(simp_m1);
		}

		if (!simp_m2.is_history_empty()) {
			current_mesh = M2;
			split_vertex(simp_m2);
		}
	
		count++;
		if (count % 100 == 0) printf("It. %d \n", count);
	}
	printf("DONE.\n");
}

void Mesh_mapping::split_vertex(Mesh_simplify& simp_m) 
{
	// ** FOR THE SAKE OF DOCUMENTATION **
	// ** SUPPOSE WE ARE VERTEX SPLITTING ON M1 **
	// ** AND OPTIMIZING ITS POSITION ON M2 **

	// Get M1 Map
	std::map<int, Inter_data>& ISM_map_1 = current_ISM(); // M1 ISM
	std::map<int, int>& ftv_map_1 = current_FTV(); // M1F TO M2V

	// Get M2 Map
	switch_mesh();
	std::map<int, Inter_data>& ISM_map_2 = current_ISM(); // M2 ISM
	std::map<int, int>& ftv_map_2 = current_FTV(); // M2F TO M1V
	switch_mesh();

	// M1 vertex split data
    HistoryEntry entry = simp_m.get_history();
	// faces that existed in M1 prior to vertex split
    std::vector<int> faces = entry.kept_faces; 
	// compute inter surface data for new vertex
    Inter_data new_v_data = ISM_iteration(entry);
		
	// get verts in M2 affected by original M1 faces
	std::vector<std::pair<int, Eigen::Vector3d>> verts_for_faces;    
	for (const auto& face : faces) {
		// check if there actually are verts in said face
		if (ftv_map_1.count(face)>0) 
		{
			int affected_v = ftv_map_1[face]; // M2 VERT
			if (ISM_map_2.count(affected_v)==0) printf("NO mapping for this vert. \n");
			Inter_data affected_v_data = ISM_map_2[affected_v]; // DATA IN M1
			Eigen::Vector3d pos = get_pos_from_inter_data(affected_v_data); // VPOS IN M1
			verts_for_faces.push_back(std::make_pair(affected_v,pos));
		}
    }

	// actually add new vert
    simp_m.construct_edge_from_history(entry);

    // NOW WE WANT TO UPDATE M2 MAP BASED ON M1 CHANGES

	// get added newly face ids
    std::vector<int> new_faces = simp_m.get_new_faces();
	// construct list of all faces in split v neighbourhood
    new_faces.insert(new_faces.end(), faces.begin(), faces.end());

    for (const auto& pair : verts_for_faces) {
        Inter_data v_data = find_enclosing_face(new_faces, pair.second); // DATA ON M1
     	ISM_map_2[pair.first] = v_data; // M2 ISM
		ftv_map_1[v_data.f_id] = pair.first; // M1F TO M2V
    }

	// get added verts id
    int new_vid = simp_m.get_split_vid();
	// populate M1 Map 
    ISM_map_1[new_vid] = new_v_data; 
	// also populate M2F TO M1V MAP
    ftv_map_2[new_v_data.f_id] = new_vid; 
}

//
// Compute map data (face, barycentric coords)
//
Inter_data Mesh_mapping::ISM_iteration(HistoryEntry &entry)
{
	// say M1, split vertex. 
	std::map<int, Inter_data>& ISM = current_ISM(); // M1 ISM

	// get faces before split, and position of new vertex
	std::vector<int> faces = entry.kept_faces;
	Eigen::Vector3d v_pos = entry.get_vertex_pos();  

	// figure out where the new vertex is relative to the original faces
	// find the face as well as its corresponding barycentric coordinates
	Inter_data v_data = find_enclosing_face(faces, v_pos); // M1 DATA 

	// given the face, get its vertices, and compute their corresponding position in M2 space
	int f_id = v_data.f_id;
	Mesh_connectivity::Face_iterator f = mesh().face_at(f_id);
	Mesh_connectivity::Half_edge_iterator he = f.half_edge();
	Mesh_connectivity::Vertex_iterator v = he.origin();
	std::vector<int> Mj_faces;
	std::vector<Eigen::Vector3d> Mjp_verts;
	switch_mesh(); // M2 
	do
	{
		Mesh_connectivity::Vertex_iterator v_curr = he.origin();
		if (ISM.count(v_curr.index())==0) printf("no mapping for this vert");
		Inter_data curr_data = ISM[v_curr.index()]; // M1 ISM -> M2 DATA
		Mj_faces.push_back(curr_data.f_id);
		Mjp_verts.push_back(get_pos_from_inter_data(curr_data));
		he = he.next();
	} while (he.origin().index() != v.index());

	// using the previously computed barycentric coordinates
	// compute the position of the new split vertex in M2 space - v'
	Eigen::Vector3d mj_vp;
	for (int i=0; i<3; i++)
	{
		mj_vp += Mjp_verts[i]*v_data.bary[i];
	}

	// for the sake of the map, now we need the face on M2, v' lives in 
	Inter_data vp_data = find_enclosing_face(Mj_faces, mj_vp);
	
	switch_mesh(); // M1 for the sake of consistency

	// return the map data
	return vp_data;
}

//
// compute position of vertex given a face and barycentric coordinates
//
Eigen::Vector3d Mesh_mapping::get_pos_from_inter_data(Inter_data data)
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

	Eigen::Vector3d v_out = Eigen::Vector3d::Zero();
	for (int i=0; i<3; i++)
	{
		v_out += verts[i] * data.bary[i];
	}
	return v_out;
}

//
// given a vertex position and a list of faces
// compute the barycentric coordinates for the vertex per face
// if we find a non-nan set of coords, return that face and corresponding coords
//
Inter_data Mesh_mapping::find_enclosing_face(std::vector<int> faces, Eigen::Vector3d v_pos)
{
	std::set<int> faces_set(faces.begin(), faces.end()); // remove any duplicates

    int f_id;
	std::vector<double> bary;

  	for (auto it = faces_set.begin(); it != faces_set.end(); ++it) {
        f_id = *it;
		bary = compute_bc(f_id, v_pos);
		if (bary[0] != -1.0) {break;}
    }
	if (bary[0] == -1.0) {printf("BAD! DIDNT FIND ENCLOSING FACE. \n");}

	return Inter_data{bary, f_id};
}

//
// Initialize Inter-surface maps given simplicial correspondence
//
void Mesh_mapping::init_ISM()
{
	for(auto &pair: mesh_map){
		current_mesh = M2;
		Inter_data data = trivial_map_data(pair.second);
		m2_ftv_map[data.f_id] = pair.first;
		ISM_M1[pair.first] = data;
		
		current_mesh = M1;
		data = trivial_map_data(pair.first);
		m1_ftv_map[data.f_id] = pair.second;
		ISM_M2[pair.second] = data;
	}
}

//
// Find the face and corresponding barycentric coordinates for a vertex
// since this is the base correspondence, coords must weight one vertex entirely
//
Inter_data Mesh_mapping::trivial_map_data(int vid)
{
	std::vector<double> bary; 

	Mesh_connectivity::Vertex_iterator vert = mesh().vertex_at(vid);
	Mesh_connectivity::Half_edge_iterator he_v = vert.half_edge();
	Mesh_connectivity::Face_iterator f_v = he_v.face();

	Mesh_connectivity::Half_edge_iterator he = f_v.half_edge();
	Mesh_connectivity::Vertex_iterator v = he.origin();
	do
	{
		Mesh_connectivity::Vertex_iterator v_curr = he.origin();
		if (v_curr.index() == vid) bary.push_back(1.0); // this is our correspondence, so weight=1
		else bary.push_back(0.0); // else weight=0
		he = he.next();
	} while (he.origin().index() != v.index());

	int f_id = f_v.index(); 

	return Inter_data{bary, f_id};
}


//
// compute the barycentric coordinates of v_pos wihin f_id
// return a default error vector if any coordinates are nan
//
std::vector<double> Mesh_mapping::compute_bc(int f_id, Eigen::Vector3d v_pos)
{
	Mesh_connectivity::Face_iterator f = mesh().face_at(f_id);
	Mesh_connectivity::Half_edge_iterator he = f.half_edge();
	Mesh_connectivity::Vertex_iterator v = he.origin();
	std::vector<double> bary;
	std::vector<double> fail;
	fail.push_back(-1.0);

	// get vertex positions on the given face
	std::vector<Eigen::Vector3d> verts;
	do
	{
		Mesh_connectivity::Vertex_iterator v_curr = he.origin();
		verts.push_back(v_curr.xyz());
		he = he.next();
	} while (he.origin().index() != v.index());

	// compute weight from the given vert to each vert on the face
	double sum = 0.0;
    for (const auto& vert : verts) 
	{
		double w = get_wij(v_pos, vert, verts);
		if (isnan(w)) {
			return fail;
		}
		else {
			bary.push_back(w);
			sum += w;
		}
	}

	// normalize
    for (auto& element : bary) {
        element /= sum;
    }

	return bary;
}

//
// compute barycentric coordinates (mean-value) between i and j
//
double Mesh_mapping::get_wij(Eigen::Vector3d i, Eigen::Vector3d j, std::vector<Eigen::Vector3d> &verts) 
{
	// Assuming verts always has size 3
    Eigen::Vector3d p, q;

    // Find the other two vertices in verts
    for (const auto& vert : verts) {
        if (vert != j) {
            if (p == Eigen::Vector3d::Zero()) {
                p = vert;
            } else {
                q = vert;
                break;
            }
        }
    }

	Eigen::Vector3d IJ = j-i;
	Eigen::Vector3d IP = p-i;
	Eigen::Vector3d IQ = q-i;

	double alpha = get_angle(IJ, IP);
	double beta = get_angle(IJ, IQ);
	double r = IJ.norm();

	return (std::tan(alpha / 2) + std::tan(beta / 2)) / r;
}

//
// update position on mesh using inter-surface map
//
void Mesh_mapping::update_positions(int _mesh) {
	current_mesh = _mesh;
	std::map<int, Inter_data>& ISM = current_ISM(); // M1 ISM

	for(auto &pair: ISM)
	{
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(pair.first); // M1 V
		switch_mesh();
		v.data().xyz = get_pos_from_inter_data(pair.second); // M2 DATA
		switch_mesh();
	}
}


double Mesh_mapping::get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) {
	return std::acos(v1.dot(v2) / (v1.norm() * v2.norm()));
}

} // end of mohe
} // end of minimesh
