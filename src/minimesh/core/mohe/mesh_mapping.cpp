#include <minimesh/core/mohe/mesh_mapping.hpp>
#include <minimesh/core/mohe/mesh_simplify.hpp>
#include <minimesh/core/util/assert.hpp>
#include <set>
#include <iostream>

namespace minimesh {
    namespace mohe {

        void logAction(std::string msg) {
            std::cout << "[ACTION] " << msg << std::endl;
        }

        void logStatus(std::string msg, bool indent = false) {
            if (indent) {
                std::cout << "\t[STATUS] " << msg << std::endl;
            } else {
                std::cout << "[STATUS] " << msg << std::endl;
            }
        }

        Mesh_mapping::VertexOptimizer Mesh_mapping::vertexOptimizer() {
            Mesh_mapping::VertexOptimizer optimizer = Mesh_mapping::VertexOptimizer(this);
            return optimizer;
        }

//
// Build Inter-surface map between two meshes
//
        void Mesh_mapping::build_mapping() {
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
            while (!simp_m1.is_history_empty() || !simp_m2.is_history_empty()) {
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

        void Mesh_mapping::split_vertex(Mesh_simplify &simp_m) {
            // ** FOR THE SAKE OF DOCUMENTATION **
            // ** SUPPOSE WE ARE VERTEX SPLITTING ON M1 **
            // ** AND OPTIMIZING ITS POSITION ON M2 **

            // Get M1 Map
            std::map<int, Inter_data> &ISM_map_1 = current_ISM(); // M1 ISM
            std::map<int, int> &ftv_map_1 = current_FTV(); // M1F TO M2V

            // Get M2 Map
            switch_mesh();
            std::map<int, Inter_data> &ISM_map_2 = current_ISM(); // M2 ISM
            std::map<int, int> &ftv_map_2 = current_FTV(); // M2F TO M1V
            switch_mesh();

            // M1 vertex split data
            HistoryEntry entry = simp_m.get_history();
            // faces that existed in M1 prior to vertex split
            std::vector<int> faces = entry.kept_faces;
            // compute inter surface data for new vertex
            Inter_data new_v_data = ISM_iteration(entry);

            // get verts in M2 affected by original M1 faces
            std::vector<std::pair<int, Eigen::Vector3d>> verts_for_faces;
            for (const auto &face: faces) {
                // check if there actually are verts in said face
                if (ftv_map_1.count(face) > 0) {
                    int affected_v = ftv_map_1[face]; // M2 VERT
                    if (ISM_map_2.count(affected_v) == 0) printf("NO mapping for this vert. \n");
                    Inter_data affected_v_data = ISM_map_2[affected_v]; // DATA IN M1
                    Eigen::Vector3d pos = get_pos_from_inter_data(affected_v_data); // VPOS IN M1
                    verts_for_faces.push_back(std::make_pair(affected_v, pos));
                }
            }

            // actually add new vert
            simp_m.construct_edge_from_history(entry);

            // NOW WE WANT TO UPDATE M2 MAP BASED ON M1 CHANGES

            // get added newly face ids
            std::vector<int> new_faces = simp_m.get_new_faces();
            // construct list of all faces in split v neighbourhood
            new_faces.insert(new_faces.end(), faces.begin(), faces.end());

            for (const auto &pair: verts_for_faces) {
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

            std::cout<< OPT_COUNT<< std::endl;

            if(OPT_COUNT > 10){
                Mesh_mapping::VertexOptimizer optimizer = vertexOptimizer();
                int sm = current_mesh == M1 ? M1 : M2;
                int tm = current_mesh == M1 ? M2 : M1;
                OutputData opt_res = optimizer.optimizeVertex(InputData{sm, tm, new_vid});

                if (opt_res.state) {
                    std::cout << "optimzied pos: " << opt_res.pos << std::endl;
                    new_v_data = find_enclosing_face(new_faces, opt_res.pos); // Inter data for optimized position
                    // populate M1 Map
                    ISM_map_1[new_vid] = new_v_data;
                }
            }

            OPT_COUNT++;
        }

//
// Compute map data (face, barycentric coords)
//
        Inter_data Mesh_mapping::ISM_iteration(HistoryEntry &entry) {
            // say M1, split vertex.
            std::map<int, Inter_data> &ISM = current_ISM(); // M1 ISM

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
            do {
                Mesh_connectivity::Vertex_iterator v_curr = he.origin();
                if (ISM.count(v_curr.index()) == 0) printf("no mapping for this vert");
                Inter_data curr_data = ISM[v_curr.index()]; // M1 ISM -> M2 DATA
                Mj_faces.push_back(curr_data.f_id);
                Mjp_verts.push_back(get_pos_from_inter_data(curr_data));
                he = he.next();
            } while (he.origin().index() != v.index());

            // using the previously computed barycentric coordinates
            // compute the position of the new split vertex in M2 space - v'
            Eigen::Vector3d mj_vp;
            for (int i = 0; i < 3; i++) {
                mj_vp += Mjp_verts[i] * v_data.bary[i];
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
        Eigen::Vector3d Mesh_mapping::get_pos_from_inter_data(Inter_data data) {
            int f_id = data.f_id;
            Mesh_connectivity::Face_iterator f = mesh().face_at(f_id);
            Mesh_connectivity::Half_edge_iterator he = f.half_edge();
            Mesh_connectivity::Vertex_iterator v = he.origin();
            std::vector<Eigen::Vector3d> verts = {};
            do {
                Mesh_connectivity::Vertex_iterator v_curr = he.origin();
                verts.push_back(v_curr.xyz());
                he = he.next();
            } while (he.origin().index() != v.index());

            Eigen::Vector3d v_out = Eigen::Vector3d::Zero();
            for (int i = 0; i < 3; i++) {
                v_out += verts[i] * data.bary[i];
            }
            return v_out;
        }

//
// given a vertex position and a list of faces
// compute the barycentric coordinates for the vertex per face
// if we find a non-nan set of coords, return that face and corresponding coords
//
        Inter_data Mesh_mapping::find_enclosing_face(std::vector<int> faces, Eigen::Vector3d v_pos) {
            std::set<int> faces_set(faces.begin(), faces.end()); // remove any duplicates

            int f_id;
            std::vector<double> bary;

            for (auto it = faces_set.begin(); it != faces_set.end(); ++it) {
                f_id = *it;
                bary = compute_bc(f_id, v_pos);
                if (bary[0] != -1.0) { break; }
            }
            if (bary[0] == -1.0) { printf("BAD! DIDNT FIND ENCLOSING FACE. \n"); }

            return Inter_data{bary, f_id};
        }

//
// Initialize Inter-surface maps given simplicial correspondence
//
        void Mesh_mapping::init_ISM() {
            for (auto &pair: mesh_map) {
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
        Inter_data Mesh_mapping::trivial_map_data(int vid) {
            std::vector<double> bary;

            Mesh_connectivity::Vertex_iterator vert = mesh().vertex_at(vid);
            Mesh_connectivity::Half_edge_iterator he_v = vert.half_edge();
            Mesh_connectivity::Face_iterator f_v = he_v.face();

            Mesh_connectivity::Half_edge_iterator he = f_v.half_edge();
            Mesh_connectivity::Vertex_iterator v = he.origin();
            do {
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
        std::vector<double> Mesh_mapping::compute_bc(int f_id, Eigen::Vector3d v_pos) {
            Mesh_connectivity::Face_iterator f = mesh().face_at(f_id);
            Mesh_connectivity::Half_edge_iterator he = f.half_edge();
            Mesh_connectivity::Vertex_iterator v = he.origin();
            std::vector<double> bary;
            std::vector<double> fail;
            fail.push_back(-1.0);

            // get vertex positions on the given face
            std::vector<Eigen::Vector3d> verts;
            do {
                Mesh_connectivity::Vertex_iterator v_curr = he.origin();
                verts.push_back(v_curr.xyz());
                he = he.next();
            } while (he.origin().index() != v.index());

            // compute weight from the given vert to each vert on the face
            double sum = 0.0;
            for (const auto &vert: verts) {
                double w = get_wij(v_pos, vert, verts);
                if (isnan(w)) {
                    return fail;
                } else {
                    bary.push_back(w);
                    sum += w;
                }
            }

            // normalize
            for (auto &element: bary) {
                element /= sum;
            }

            return bary;
        }

//
// compute barycentric coordinates (mean-value) between i and j
//
        double Mesh_mapping::get_wij(Eigen::Vector3d i, Eigen::Vector3d j, std::vector<Eigen::Vector3d> &verts) {
            // Assuming verts always has size 3
            Eigen::Vector3d p, q;

            // Find the other two vertices in verts
            for (const auto &vert: verts) {
                if (vert != j) {
                    if (p == Eigen::Vector3d::Zero()) {
                        p = vert;
                    } else {
                        q = vert;
                        break;
                    }
                }
            }

            Eigen::Vector3d IJ = j - i;
            Eigen::Vector3d IP = p - i;
            Eigen::Vector3d IQ = q - i;

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
            std::map<int, Inter_data> &ISM = current_ISM(); // M1 ISM

            for (auto &pair: ISM) {
                Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(pair.first); // M1 V
                switch_mesh();
                v.data().xyz = get_pos_from_inter_data(pair.second); // M2 DATA
                switch_mesh();
            }
        }


        double Mesh_mapping::get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) {
            return std::acos(v1.dot(v2) / (v1.norm() * v2.norm()));
        }


        /*
         * ================================================================
         * VERTEX OPTIMIZER
         * =================================================================
         *
         * */

        void Mesh_mapping::VertexOptimizer::prime_optimizer() {

            std::cout << "Priming Optimizer" << std::endl;

            Mesh_connectivity::Vertex_ring_iterator source_ring = ring_at_source();
            source_vertex = ism->mesh().vertex_at(_opt_data.source_v_id);

            // collect the ring data for the source vertex
            logAction("Collecting Source Ring");
            do {
                v_map[source_ring.half_edge().origin().index()] = source_ring.half_edge().origin().xyz();
                std::cout <<source_ring.half_edge().origin().xyz() <<std::endl;
            } while (source_ring.advance());
            logStatus("successful", true);

            logAction("Collecting Image Ring");
            // collect the ring data for the target vertex

            Mesh_connectivity::Vertex_ring_iterator source_ring2 = ring_at_source();
            do {
                v_bar_map[source_ring2.half_edge().origin().index()] = get_image_position(
                        source_ring2.half_edge().origin().index());
            } while (source_ring2.advance());
            logStatus("Complete", true);

            //  compute the flattened positions
            int s = static_cast<int>(v_bar_map.size());
            std::vector<Eigen::Vector3d> vectors_v_bar = {};

            // collect the vectors for vw_bar
            std::vector<double> lengths = {};
            for (const auto &pair: v_bar_map) {
                vectors_v_bar.push_back(pair.second - source_vertex.xyz());
                lengths.push_back((pair.second - source_vertex.xyz()).norm());
            }

            double baseSectorSize = 2 * numbers::pi / s;
            // compute the angles between them, collect the ratios
            std::vector<double> v_bar_angles = {};
            double v_bar_angle_sum = 0;
            for (int i = 0; i < s; ++i) {
                int j = (i + 1) % s;
                double _angle = acos((vectors_v_bar[i].dot(vectors_v_bar[j])) /
                                     (vectors_v_bar[i].norm() * vectors_v_bar[j].norm()));
                v_bar_angles.push_back(_angle);
                v_bar_angle_sum += _angle;
            }

            // collect the ring v_hat_positions,
            std::vector<Eigen::Vector3d> ring_v_hat_positions = {};
            for (int i = 0; i < v_bar_angles.size(); i++) {
                double th = v_bar_angles[i] * (360 / v_bar_angle_sum);
                double x = lengths[i] * cos(th);
                double y = lengths[i] * sin(th);
                ring_v_hat_positions.push_back(Eigen::Vector3d(x, y, 0.0));
            }

            // assign the collected positions to their map
            int i = 0;
            for (const auto &pair: v_bar_map) {
                v_hat_map[pair.first] = ring_v_hat_positions[i];
                i++;
            }

            logStatus("priming complete");
        }

        void Mesh_mapping::VertexOptimizer::relax_neighbourhood() {
            logAction("Relaxing neighbourhood");
            // create matrix system and solve for free boundary parameterization
            // adapt code from free_param_to_here

            // define the holders for the algorithm below

            // list of boundary vertices
            std::vector<int> boundary_ids;

            // map from MF1/MF2 column index to vertex id
            std::map<int, int> free;
            std::map<int, int> free_rev; //reverse

            // map of new vertex positions
            std::map<int, Eigen::Vector3d> new_positions;

            // Index of boundary vertices
            int num_boundary;

            // vid -> mat id
            std::map<int, int> interior;

            // mat id -> vid
            std::map<int, int> interior_rev;

            // map of new vertex positions

            Eigen::MatrixXd Ubar;
            Eigen::MatrixXd Vbar;
            Eigen::SparseMatrix<double> A;
            std::vector<Eigen::Triplet<double>> A_elem;

            // locate points in v_bar that exists in the boundary formed by v_hat
            // solve the parameterization least squares system

            // create list of active boundary points
            std::vector<Eigen::Vector3d> boundaryPoints = {};
            for (const auto &pair: v_hat_map) {
                boundaryPoints.push_back(pair.second);
                boundary_ids.push_back(pair.first);
            }

            // collect the selected points
            for (const auto &pair: v_bar_map) {
                if (isPointInsideBoundary(boundaryPoints, pair.second)) {
                    chosenPoints.push_back(pair.second);
                    chosenPointsIds.push_back(pair.first);
                }
            }

            // faces belonging to points inside the boundary
            // map of the half edge and its corresponding face
            Mesh_connectivity::Vertex_ring_iterator source_ring3 = ring_at_source();

            do {
                int idx = source_ring3.half_edge().origin().index();
                ArraySearchResults res = arrayContains(chosenPointsIds, idx);
                if (res.found) {
                    chosenFaces[source_ring3.half_edge().index()] = source_ring3.half_edge().face().index();
                }
            } while (source_ring3.advance());

            // construct the mass sparse matrix

            /*
             * Prime the components
             * 1. interior
             * 2. boundary values
             * 3. number of boundary vertices
             * */
            {
                for (int i = 0; i < static_cast<int>(chosenPointsIds.size()); ++i) {
                    interior[chosenPointsIds[i]] = i;
                    interior_rev[i] = chosenPointsIds[i];
                }

                for (const auto &pair: v_hat_map) {
                    new_positions[pair.first] = pair.second;
                }

            }

            /* perform math for all matrix values*/

            logStatus("\tconstructing free boundary param matrix", true);
            {
                int n = static_cast<int>(chosenPoints.size());
                A.resize(n, n);
                A.setZero();
                Ubar = Eigen::MatrixXd::Zero(n, 1);
                Vbar = Eigen::MatrixXd::Zero(n, 1);
                for (const auto &pair: interior) {
                    int vid = pair.first;
                    int i = pair.second;
                    // compute A_i
                    {
                        Mesh_connectivity::Vertex_ring_iterator ring = get_ring_from_target_mesh(vid);

                        do {
                            Mesh_connectivity::Vertex_iterator v_j = ring.half_edge().origin();

                            // check adjacent vertex is NOT on the boundary
                            ArraySearchResults res = arrayContains(boundary_ids, vid);
                            if (!res.found) {
                                int jid = v_j.index();
                                int j = interior[jid];

                                // A_{ij} is inverse normalized weight for ij
                                A_elem.push_back(Eigen::Triplet<double>(i, j, -lambda_ij(vid, jid)));
                            }
                        } while (ring.advance());

                    }

                    // compute_UVbar_i
                    {
                        Mesh_connectivity::Vertex_ring_iterator ring = get_ring_from_target_mesh(vid);
                        double u_sum = 0.0;
                        double v_sum = 0.0;
                        do // *points to*
                        {
                            Mesh_connectivity::Vertex_iterator v_j = ring.half_edge().origin();

                            // check adjacent vertex IS on the boundary
                            if (v_j.is_boundary()) {
                                int jid = v_j.index();

                                // normalized weight for ij
                                double k_ij = lambda_ij(vid, jid);

                                // new unit circle positions for j
                                double u_j = new_positions[jid][0];
                                double v_j = new_positions[jid][1];

                                double u_ij = k_ij * u_j;
                                double v_ij = k_ij * v_j;

                                u_sum += u_ij;
                                v_sum += v_ij;
                            }
                        } while (ring.advance());

                        Ubar(i) = u_sum;
                        Vbar(i) = v_sum;
                    }
                }
                A.setFromTriplets(A_elem.begin(), A_elem.end());
            }


            logAction("Solving linear system");
            /* solve the linear system */
            {
                // Use LU since A is not symmetric!
                Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
                solver.analyzePattern(A);
                solver.factorize(A);
                Eigen::MatrixXd U = solver.solve(Ubar);
                Eigen::MatrixXd V = solver.solve(Vbar);

                // populate vertex position map with solved UVs
                for (int i = 0; i < U.rows(); ++i) {
                    new_positions[interior_rev[i]] = Eigen::Vector3d(U(i), V(i), 0.0);
                }
            }

            /* collect the results of the flattened section */
            {
                relaxed_v_hat_map = new_positions;
            }


            logStatus("Relaxation Complete");
        }

        void Mesh_mapping::VertexOptimizer::perform_line_with_distortion_calculations() {

            logAction("Performing line search");

            // vertex whose position we want to optimize
            Eigen::Vector3d v = source_vertex.xyz();
            Eigen::Vector3d optimized_v = source_vertex.xyz();

            // convergence tolerance for the line search
            double tolerance = 1.0e-5;

            // original distortion
            double o_distortion = 0.0;
            double min_distortion = 10.0; // define an unusually high distortion

            // distortions computed per line search
            std::vector<double> distortions;

            //
            // begin optimization
            //

            // create list of v_hat_kernel
            std::vector<Eigen::Vector3d> v_hat_kernel;
            for (auto &pair: v_hat_map) {
                v_hat_kernel.push_back(pair.second);
            }

            double c_error = 1.0;
            int n_iterations = 200;

            logStatus("starting line search iterations", true);
            while (c_error > tolerance) {
                // perform line search and compute distortion

                /*
                 * Perform line search optimization
                 * */
                {
                    // chose a random direction of v inside v_hat_kernel
                    int idx_1 = getRandomValue(static_cast<int>(v_hat_kernel.size()));
                    int idx_2 = getRandomValue(static_cast<int>(v_hat_kernel.size()));
                    Eigen::Vector3d dir = v_hat_kernel[idx_1] - v_hat_kernel[idx_2];
                    double t = random_zero_to_one();
                    optimized_v = optimized_v + (dir * t);
                    double pm_error = compute_pm_error_over_target(optimized_v);
                    c_error = pm_error - c_error;
                }

                /*
                 * Compute distortion metric
                 * */
                {
                    double stretch_distortion = 0.0;
                    double area_m1 = 0.0;
                    double area_m2 = 0.0;
                    double L = 0.0;
                    double G = 0.0;

                    // for every triangle, compute the left and right distortion metric, and combine
                    //   get corresponding vertex positions from map

                    /*
                     * Compute area of region M1
                     * */
                    {
                        Mesh_connectivity::Vertex_ring_iterator source_ring = ring_at_source();

                        do {
                            Mesh_connectivity::Face_iterator face = source_ring.half_edge().face();
                            Mesh_connectivity::Vertex_iterator v1 = face.half_edge().origin();
                            Mesh_connectivity::Vertex_iterator v2 = face.half_edge().dest();
                            Mesh_connectivity::Vertex_iterator v3 = face.half_edge().next().dest();
                            Eigen::Vector3d q1 = v1.xyz();
                            Eigen::Vector3d q2 = v2.xyz();
                            Eigen::Vector3d q3 = v3.xyz();

                            area_m1 = area_m1 + (q1.dot(q2) + q1.dot(q3) + q2.dot(q3));

                        } while (source_ring.advance());
                    }

                    /*
                     * Compute the area or region M2 (its image)
                     * */
                    {

                        Mesh_connectivity::Vertex_ring_iterator source_ring = ring_at_source();
                        do {
                            Mesh_connectivity::Face_iterator face = source_ring.half_edge().face();
                            Mesh_connectivity::Vertex_iterator v1 = face.half_edge().origin();
                            Mesh_connectivity::Vertex_iterator v2 = face.half_edge().dest();
                            Mesh_connectivity::Vertex_iterator v3 = face.half_edge().next().dest();
                            Eigen::Vector3d q1 = get_image_position(v1.index());
                            Eigen::Vector3d q2 = get_image_position(v2.index());
                            Eigen::Vector3d q3 = get_image_position(v3.index());

                            area_m2 = area_m2 + (q1.dot(q2) + q1.dot(q3) + q2.dot(q3));

                        } while (source_ring.advance());
                    }

                    /*
                     * Compute L and G
                     * */

                    {
                        /* get the target coordinate in via the map, and used the flattened u,v coordinates */
                        Mesh_connectivity::Face_iterator face;
                            face = ism->mesh().face_at(ism->ISM_M1[_opt_data.source_v_id].f_id);

                        // define points on the face
                        Mesh_connectivity::Vertex_iterator v1 = face.half_edge().origin();
                        Mesh_connectivity::Vertex_iterator v2 = face.half_edge().dest();
                        Mesh_connectivity::Vertex_iterator v3 = face.half_edge().next().dest();
                        Eigen::Vector3d q1 = get_image_position(v1.index());
                        Eigen::Vector3d q2 = get_image_position(v2.index());
                        Eigen::Vector3d q3 = get_image_position(v3.index());

                        // compute L and G

                        Eigen::Vector3d p1 = relaxed_v_hat_map[v1.index()];
                        Eigen::Vector3d p2 = relaxed_v_hat_map[v2.index()];
                        Eigen::Vector3d p3 = relaxed_v_hat_map[v3.index()];

                        double area = q1.dot(q2) + q1.dot(q3) + q2.dot(q3);

                        Eigen::Vector3d S_s =
                                (q1 * (p2.y() - p3.y()) + q2 * (p3.y() - p1.y()) + q3 * (p1.y() - p2.y())) /
                                (2 * area);
                        Eigen::Vector3d S_t =
                                (q1 * (p3.x() - p2.x()) + q2 * (p1.x() - p3.x()) + q3 * (p2.x() - p1.x())) /
                                (2 * area);

                        double a = S_s.dot(S_s);
                        double b = S_s.dot(S_t);
                        double c = S_t.dot(S_t);

                        L = sqrt(0.5 * ((a + c) + sqrt(pow(a - c, 2) + (4 * pow(b, 2)))));
                        G = sqrt(0.5 * ((a + c) - sqrt(pow(a - c, 2) + (4 * pow(b, 2)))));

                    }


                    /*
                     * Compute Left and right triangle distortions and sum them
                     * */

                    {
                        double left_sum = 0;
                        double right_sum = 0;

                        Mesh_connectivity::Vertex_ring_iterator source_ring = ring_at_source();

                        do {
                            Mesh_connectivity::Face_iterator face = source_ring.half_edge().face();
                            Mesh_connectivity::Vertex_iterator v1 = face.half_edge().origin();
                            Mesh_connectivity::Vertex_iterator v2 = face.half_edge().dest();
                            Mesh_connectivity::Vertex_iterator v3 = face.half_edge().next().dest();
                            Eigen::Vector3d q1 = get_image_position(v1.index());
                            Eigen::Vector3d q2 = get_image_position(v2.index());
                            Eigen::Vector3d q3 = get_image_position(v3.index());

                            double area = q1.dot(q2) + q1.dot(q3) + q2.dot(q3);
                            double _dis = area * (area_m1 / pow(area_m2, 2)) * (pow(G, 2) + pow(L, 2));
                            right_sum = right_sum + _dis;
                        } while (source_ring.advance());

                        Mesh_connectivity::Vertex_ring_iterator source_ring2 = ring_at_source();
                        do {
                            Mesh_connectivity::Face_iterator face = source_ring2.half_edge().face();
                            Mesh_connectivity::Vertex_iterator v1 = face.half_edge().origin();
                            Mesh_connectivity::Vertex_iterator v2 = face.half_edge().dest();
                            Mesh_connectivity::Vertex_iterator v3 = face.half_edge().next().dest();
                            Eigen::Vector3d q1 = v1.xyz();
                            Eigen::Vector3d q2 = v2.xyz();
                            Eigen::Vector3d q3 = v3.xyz();

                            double area = q1.dot(q2) + q1.dot(q3) + q2.dot(q3);
                            double _dis = area * (area_m2 / pow(area_m1, 2)) * ((1 / pow(G, 2)) + (1 / pow(L, 2)));
                            left_sum = left_sum + _dis;
                        } while (source_ring2.advance());

                        stretch_distortion = left_sum + right_sum;
                    }

                    // get the smallest distortion
                    if (stretch_distortion < min_distortion) min_distortion = stretch_distortion;
                }
                n_iterations--;
                if (n_iterations == 0) break;
            }

            // for now, do not check previous and current distortions
            /*if (min_distortion >= o_distortion) {
                _res.state = false;
                _res.pos = v;
            } else {
                _res.state = true;
                _res.pos = optimized_v;
            }*/

            _res.state = true;
            _res.pos = optimized_v;
        }

        OutputData Mesh_mapping::VertexOptimizer::optimizeVertex(const minimesh::mohe::InputData &data) {
            _opt_data = data;

            // prime the optimizer
            prime_optimizer();

            // relax points to 2d space
            relax_neighbourhood();

            // line search with distortion
            perform_line_with_distortion_calculations();

            return _res;
        }

        Eigen::Vector3d Mesh_mapping::VertexOptimizer::get_image_position(int v_id) {
            return ism->get_pos_from_inter_data(ism->current_ISM()[v_id]);
        }

        double Mesh_mapping::VertexOptimizer::compute_pm_error_over_target(Eigen::Vector3d &v) {
                return compute_pm_error_over_target(v, ism->mesh());
        }

        double
        Mesh_mapping::VertexOptimizer::compute_pm_error_over_target(Eigen::Vector3d &v, Mesh_connectivity &mesh) {
            double error = 0.0;
            for (auto &pair: v_hat_map) {
                error = error + error_over_neighbourhood(pair.first, pair.second, mesh);
            }
            return error;
        }

        double Mesh_mapping::VertexOptimizer::error_over_neighbourhood(int w_index, Eigen::Vector3d &v,
                                                                       Mesh_connectivity &mesh) {
            double lambda = 0.5;
            double psi = 1 / static_cast<double>(relaxed_v_hat_map.size());
            double deviation_error = 1.0;
            double stretch_error = 0.0;
            double sum_area = 0.0;
            // compute area for the flattened section
            for (auto &pair: chosenFaces) {
                Mesh_connectivity::Face_iterator face = get_face_from_target_mesh(pair.second);
                Mesh_connectivity::Vertex_iterator v1 = face.half_edge().origin();
                Mesh_connectivity::Vertex_iterator v2 = face.half_edge().dest();
                Mesh_connectivity::Vertex_iterator v3 = face.half_edge().next().dest();
                Eigen::Vector3d q1 = v1.xyz();
                Eigen::Vector3d q2 = v2.xyz();
                Eigen::Vector3d q3 = v3.xyz();
                Eigen::Vector3d p1 = relaxed_v_hat_map[v1.index()];
                Eigen::Vector3d p2 = relaxed_v_hat_map[v2.index()];
                Eigen::Vector3d p3 = relaxed_v_hat_map[v3.index()];

                double area = (((p2(0) - p1(0)) * (p3(1) - p1(1))) - ((p3(0) - p1(0)) * (p2(1) - p1(1)))) / 2;

                Eigen::Vector3d S_s =
                        (q1 * (p2.y() - p3.y()) + q2 * (p3.y() - p1.y()) + q3 * (p1.y() - p2.y())) / (2 * area);
                Eigen::Vector3d S_t =
                        (q1 * (p3.x() - p2.x()) + q2 * (p1.x() - p3.x()) + q3 * (p2.x() - p1.x())) / (2 * area);

                double a = S_s.dot(S_s);
                double b = S_s.dot(S_t);
                double c = S_t.dot(S_t);

                double se = sqrt(0.5 * (a + c));

                sum_area = sum_area + area;
                stretch_error = stretch_error + se;
            }

            return psi *
                   (lambda * (stretch_error) + (1 - lambda) * (deviation_error * deviation_error) / sum_area);
        }

        bool
        Mesh_mapping::VertexOptimizer::isPointInsideBoundary(const std::vector<Eigen::Vector3d> &boundary,
                                                             const Eigen::Vector3d &point) {
            int intersectionCount = 0;

            // Iterate over each edge of the boundary
            for (size_t i = 0; i < boundary.size(); ++i) {
                Eigen::Vector3d p1 = boundary[i];
                Eigen::Vector3d p2 = boundary[(i + 1) % boundary.size()];

                // Check if the line segment intersects with the horizontal line from the point
                if ((p1.y() < point.y() && p2.y() >= point.y()) ||
                    (p2.y() < point.y() && p1.y() >= point.y())) {
                    // Calculate the x-coordinate of the intersection point
                    double intersectionX = (p1.x() +
                                            (point.y() - p1.y()) / (p2.y() - p1.y()) * (p2.x() - p1.x()));

                    // Check if the intersection point is to the right of the point
                    if (point.x() < intersectionX) {
                        intersectionCount++;
                    }
                }

                // Check if the point lies on a vertex of the boundary
                if ((p1.y() == point.y() && p1.x() == point.x()) ||
                    (p2.y() == point.y() && p2.x() == point.x())) {
                    return false; // Point lies on the boundary, consider it outside
                }
            }

            // If the number of intersections is odd, the point is inside the boundary
            return (intersectionCount % 2 == 1);
        }

        Mesh_connectivity::Vertex_ring_iterator Mesh_mapping::VertexOptimizer::get_ring_from_target_mesh(int v_id) {
            return ism->mesh().vertex_ring_at(v_id);
        }

        Mesh_connectivity::Face_iterator Mesh_mapping::VertexOptimizer::get_face_from_target_mesh(int f_id) {
                return ism->mesh().face_at(f_id);
        }

        double Mesh_mapping::VertexOptimizer::lambda_ij(int i, int j) {
                return lambda_ij_target(i, j, ism->mesh());
        }

        double Mesh_mapping::VertexOptimizer::get_wik_target(int he_index, Mesh_connectivity &mesh) {
            // points to i from j
            Mesh_connectivity::Half_edge_iterator he = mesh.half_edge_at(he_index);

            Mesh_connectivity::Vertex_iterator I = he.dest();
            Mesh_connectivity::Vertex_iterator J = he.origin();
            Mesh_connectivity::Vertex_iterator P = he.twin().next().dest();
            Mesh_connectivity::Vertex_iterator Q = he.next().dest();

            Eigen::Vector3d I_xyz = I.xyz();
            Eigen::Vector3d J_xyz = J.xyz();
            Eigen::Vector3d P_xyz = P.xyz();
            Eigen::Vector3d Q_xyz = Q.xyz();

            Eigen::Vector3d IJ = J_xyz - I_xyz;
            Eigen::Vector3d IP = P_xyz - I_xyz;
            Eigen::Vector3d IQ = Q_xyz - I_xyz;

            double alpha = get_angle(IJ, IP);
            double beta = get_angle(IJ, IQ);
            double r = IJ.norm();

            return (std::tan(alpha / 2) + std::tan(beta / 2)) / r;

        }

        double Mesh_mapping::VertexOptimizer::lambda_ij_target(int i, int j, Mesh_connectivity &mesh) {

            Mesh_connectivity::Vertex_iterator v_i = mesh.vertex_at(i);
            Mesh_connectivity::Vertex_ring_iterator ring = mesh.vertex_ring_at(i);
            double w_sum = 0.0;
            double w_ij = 0.0;

            // compute MVW for all adjacent vertices
            do // *points to*
            {
                int k = ring.half_edge().origin().index();
                double w_ik = get_wik_target(ring.half_edge().index(), mesh);

                // k is target adajcent vertex j
                if (k == j) {
                    w_ij = w_ik;
                }

                w_sum += w_ik;

            } while (ring.advance());

            return w_ij / w_sum;
        }

        Mesh_connectivity::Vertex_ring_iterator Mesh_mapping::VertexOptimizer::ring_at_source() {
            return ism->mesh().vertex_ring_at(_opt_data.source_v_id);
        }


    } // end of mohe
} // end of minimesh
