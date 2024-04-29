#include <minimesh/core/mohe/mesh_mapping.hpp>
#include <minimesh/core/mohe/mesh_simplify.hpp>
#include <minimesh/core/util/assert.hpp>
#include <set>

namespace minimesh
{
namespace mohe
{


void Mesh_mapping::build_mapping()
{
	// inter surface mapping

    Mesh_simplify simp_m1(_m1);
	simp_m1.simplify_to_target(4);
    Mesh_simplify simp_m2(_m2);
	simp_m2.simplify_to_target(4);

	// TODO NEED A MAPPING BETWEEN M1 AND M2

	while (!simp_m1.is_history_empty() || !simp_m2.is_history_empty()) 
	{
		if (!simp_m1.is_history_empty()) {
			current_mesh = M1;

			HistoryEntry entry = simp_m1.get_history();
			std::vector<int> faces = entry.kept_faces;
			Inter_data new_v_data = ISM_iteration(entry);
			simp_m1.construct_edge_from_history(entry);
			int new_vid = simp_m1.get_split_vid();
			ISM_M1[new_vid] = new_v_data;
			m2_ftv_map[new_v_data.f_id] = new_vid;

			// here we update any mappings from m2 that concern the split neighbourhood

			std::vector<int> new_faces = simp_m1.get_new_faces();

			// Combine faces and new_faces using insert function
			new_faces.insert(new_faces.end(), faces.begin(), faces.end());

			std::vector<int> verts_for_faces;
			for (const auto& face : faces) 
			{
				verts_for_faces.push_back(m1_ftv_map[face]);
			}

			switch_mesh();

			for (const auto& v_id : verts_for_faces) 
			{
				Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(v_id);
				Inter_data v_data = find_enclosing_face(new_faces, v.xyz());
				ISM_M2[v_id] = v_data;
			}
			
			switch_mesh();

		}

		if (!simp_m2.is_history_empty()) {
			current_mesh = M2;

			HistoryEntry entry = simp_m2.get_history();
			std::vector<int> faces = entry.kept_faces;
			Inter_data new_v_data = ISM_iteration(entry);
			simp_m2.construct_edge_from_history(entry);
			int new_vid = simp_m2.get_split_vid();
			ISM_M2[new_vid] = new_v_data;
			m2_ftv_map[new_v_data.f_id] = new_vid;

			// here we update any mappings from m2 that concern the split neighbourhood

			std::vector<int> new_faces = simp_m2.get_new_faces();

			// Combine faces and new_faces using insert function
			new_faces.insert(new_faces.end(), faces.begin(), faces.end());

			std::vector<int> verts_for_faces;
			for (const auto& face : faces) 
			{
				verts_for_faces.push_back(m2_ftv_map[face]);
			}

			switch_mesh();

			for (const auto& v_id : verts_for_faces) 
			{
				Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(v_id);
				Inter_data v_data = find_enclosing_face(new_faces, v.xyz());
				ISM_M1[v_id] = v_data;
			}
			
			switch_mesh();
		}
	}
}


Inter_data Mesh_mapping::ISM_iteration(HistoryEntry &entry)
{
	// say M1, split vertex. 
	std::map<int, Inter_data> & ISM = ISM_M1;
	if (current_mesh != M1) ISM = ISM_M2;

	// 2. get keptFaces and removedVertex from history Entry
	std::vector<int> faces = entry.kept_faces;
	Eigen::Vector3d v_pos = entry.get_vertex_pos();  

	// 3. compute bary coords on keptFaces vs. removedVertex to see where removedVertex lives
	Inter_data v_data = find_enclosing_face(faces, v_pos);

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
		Inter_data curr_data = ISM[v_curr.index()];
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
	Inter_data vp_data = find_enclosing_face(Mj_faces, mj_vp);

	// 8. add v to the map m1->m2' or v1 -> F, C
	// ISM[v_id] = vp_data;
	return vp_data;
}

//
// compute position of vertice using barycentric coordinates given by inter data
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

	Eigen::Vector3d v_out;
	for (int i=0; i<3; i++)
	{
		v_out += verts[i]*data.bary[i];
	}
	return v_out;
}

//
// compute the barycentric coordinates of v_pos wihin each face in faces
// return the face and coordinates of the face with non-negative bary. coords
//
Inter_data Mesh_mapping::find_enclosing_face(std::vector<int> faces, Eigen::Vector3d v_pos)
{
	std::set<int> faces_set(faces.begin(), faces.end());

    // Compute barycentric coordinates for each element in the set
    int f_id;
	std::vector<double> bary;

  	for (auto it = faces_set.begin(); it != faces_set.end(); ++it) {
        f_id = *it;
		bary = compute_bc(f_id, v_pos);
		if (bary[0] != -1.0) {break;}
    }

	return Inter_data{bary, f_id};
}


void Mesh_mapping::init_ISM()
{
	for(auto &pair: mesh_map){
		Inter_data data = trivial_map_data(pair.second);
		m1_ftv_map[data.f_id] = pair.first;
		ISM_M1[pair.first] = data;
		
		data = trivial_map_data(pair.first);
		m2_ftv_map[data.f_id] = pair.second;
		ISM_M2[pair.second] = data;
	}
}

Inter_data Mesh_mapping::trivial_map_data(int vid)
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

	return Inter_data{bary, f_id};
}


//
// compute the barycentric coordinates of v_pos wihin f_id
//
std::vector<double> Mesh_mapping::compute_bc(int f_id, Eigen::Vector3d v_pos)
{
	Mesh_connectivity::Face_iterator f = mesh().face_at(f_id);
	Mesh_connectivity::Half_edge_iterator he = f.half_edge();
	Mesh_connectivity::Vertex_iterator v = he.origin();
	std::vector<double> bary;

	std::vector<double> fail;
	fail.push_back(-1.0);

	std::vector<Eigen::Vector3d> verts;
	do
	{
		Mesh_connectivity::Vertex_iterator v_curr = he.origin();
		verts.push_back(v_curr.xyz());
		he = he.next();
	} while (he.origin().index() != v.index());

    for (const auto& vert : verts) 
	{
		double w = get_wij(v_pos, vert, verts);
		if (w < 0.0) {
			return fail;
		}
		else {
			bary.push_back(w);
		}
	}

	return bary;
}

// barycentric coordinates (mean-value)
// the weight of vertex j for i
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



} // end of mohe
} // end of minimesh
