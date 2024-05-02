////
//// Created by ndehn on 2024-04-28.
////
//
//#include <minimesh/core/mohe/mesh_connectivity.hpp>
//#include "minimesh/core/util/numbers.hpp"
//#include <algorithm>
//#include <vector>
//#include <Eigen/Dense>
//#include <Eigen/SparseLU>
//#include <string>
//#include <random>
//#include <unordered_set>
//#include <ctime>
//#include <cmath>
//
//namespace minimesh {
//    namespace mohe {
//
//        /*
//         * Structure of data for calling the optimization function in the optimizer class
//         * Source mesh is the mesh we are coming from, with default value "m1"
//         * Target mesh it the mesh we are going to, with default value "m2"
//         * v_id is the vertex whose position we want to optimize from m2 to m1
//         *
//         * */
//        struct InputData {
//            std::string source_mesh = "m2";
//            std::string target_mesh = "m1";
//            int source_v_id;
//            int target_v_id;
//        };
//
//        struct OutputData {
//            bool state = false;
//            Eigen::Vector3d pos = Eigen::Vector3d(0, 0, 0);
//        };
//
//        struct ArraySearchResults {
//            bool found = false;
//            int index = -1;
//        };
//
//        struct ComplexNumber {
//            double real;
//            double imaginary;
//        };
//
//        struct IsPinnedVertex {
//            bool pinned;
//            int index;
//        };
//
//        /*
//         * Trivial Class representing InterSurface Mapping
//         * */
//        class ISM {
//        public:
//            ISM(Mesh_connectivity &mesh1, Mesh_connectivity &mesh2) : _m1(mesh1), _m2(mesh2) {}
//
//            // return M2
//            Mesh_connectivity &m1() { return _m1; }
//
//            // return M1
//            Mesh_connectivity &m2() { return _m2; }
//
//        private:
//            Mesh_connectivity &_m1;
//            Mesh_connectivity &_m2;
//        };
//
//
//        class VertexOptimizer {
//        public:
//
//            /* default constructor */
//            VertexOptimizer() = default;
//
//            /* Get the optimization results*/
//            OutputData res() { return _res; }
//
//            /*
//             * Gets the ring for the target mesh, computes the flattened positions in 2d space and saves them their
//             * required mapped regions.
//             * */
//            void prime_optimizer() {
//
//                if (_opt_data.target_mesh == "m1") {
//                    source_ring = ism->m1().vertex_ring_at(_opt_data.source_v_id);
//                    source_vertex = ism->m1().vertex_at(_opt_data.source_v_id);
//                } else {
//                    target_ring = ism->m2().vertex_ring_at(_opt_data.target_v_id);
//                    target_vertex = ism->m2().vertex_at(_opt_data.target_v_id);
//                }
//
//                if (_opt_data.target_mesh == "m2") {
//                    source_ring = ism->m2().vertex_ring_at(_opt_data.source_v_id);
//                    source_vertex = ism->m2().vertex_at(_opt_data.source_v_id);
//                } else {
//                    target_ring = ism->m1().vertex_ring_at(_opt_data.target_v_id);
//                    target_vertex = ism->m1().vertex_at(_opt_data.target_v_id);
//                }
//
//                // collect the ring data for the source vertex
//                do {
//                    v_map[source_ring.half_edge().origin().index()] = source_ring.half_edge().origin().xyz();
//                } while (source_ring.advance());
//
//                // collect the ring data for the target vertex
//                do {
//                    v_bar_map[target_ring.half_edge().origin().index()] = target_ring.half_edge().origin().xyz();
//                } while (target_ring.advance());
//
//                //  compute the flattened positions
//                int s = static_cast<int>(v_bar_map.size());
//                std::vector<Eigen::Vector3d> vectors_v_bar = {};
//
//                // collect the vectors for vw_bar
//                std::vector<double> lengths = {};
//                for (const auto &pair: v_bar_map) {
//                    vectors_v_bar.push_back(pair.second - target_vertex.xyz());
//                    lengths.push_back((pair.second - target_vertex.xyz()).norm());
//                }
//
//                double baseSectorSize = 2 * numbers::pi / s;
//                // compute the angles between them, collect the ratios
//                std::vector<double> v_bar_angles = {};
//                double v_bar_angle_sum = 0;
//                for (int i = 0; i < s; ++i) {
//                    int j = (i + 1) % s;
//                    double _angle = acos((vectors_v_bar[i].dot(vectors_v_bar[j])) /
//                                         (vectors_v_bar[i].norm() * vectors_v_bar[j].norm()));
//                    v_bar_angles.push_back(_angle);
//                    v_bar_angle_sum += _angle;
//                }
//
//                // collect the ring v_hat_positions,
//                std::vector<Eigen::Vector3d> ring_v_hat_positions = {};
//                for (int i = 0; i < v_bar_angles.size(); i++) {
//                    double th = v_bar_angles[i] * (360 / v_bar_angle_sum);
//                    double x = lengths[i] * cos(th);
//                    double y = lengths[i] * sin(th);
//                    ring_v_hat_positions.push_back(Eigen::Vector3d(x, y, 0.0));
//                }
//
//                // assign the collected positions to their map
//                int i = 0;
//                for (const auto &pair: v_bar_map) {
//                    v_hat_map[pair.first] = ring_v_hat_positions[i];
//                    i++;
//                }
//
//            }
//
//            /*
//             * Mean value parameterization to relax points in N(v_bar) into N(v_hat)
//             * */
//            void relax_neighbourhood() {
//
//                // create matrix system and solve for free boundary parameterization
//                // adapt code from free_param_to_here
//
//                // define the holders for the algorithm below
//
//                // list of boundary vertices
//                std::vector<int> boundary_ids;
//
//                // map from MF1/MF2 column index to vertex id
//                std::map<int, int> free;
//                std::map<int, int> free_rev; //reverse
//
//                // map of new vertex positions
//                std::map<int, Eigen::Vector3d> new_positions;
//
//                // Index of boundary vertices
//                int num_boundary;
//
//                // vid -> mat id
//                std::map<int, int> interior;
//
//                // mat id -> vid
//                std::map<int, int> interior_rev;
//
//                // map of new vertex positions
//
//                Eigen::MatrixXd Ubar;
//                Eigen::MatrixXd Vbar;
//                Eigen::SparseMatrix<double> A;
//                std::vector<Eigen::Triplet<double>> A_elem;
//
//                // locate points in v_bar that exists in the boundary formed by v_hat
//                // solve the parameterization least squares system
//
//                // create list of active boundary points
//                std::vector<Eigen::Vector3d> boundaryPoints = {};
//                for (const auto &pair: v_hat_map) {
//                    boundaryPoints.push_back(pair.second);
//                    boundary_ids.push_back(pair.first);
//                }
//
//                // collect the selected points
//                for (const auto &pair: v_bar_map) {
//                    if (isPointInsideBoundary(boundaryPoints, pair.second)) {
//                        chosenPoints.push_back(pair.second);
//                        chosenPointsIds.push_back(pair.first);
//                    }
//                }
//
//                // faces belonging to points inside the boundary
//                // map of the half edge and its corresponding face
//                do {
//                    int idx = target_ring.half_edge().origin().index();
//                    ArraySearchResults res = arrayContains(chosenPointsIds, idx);
//                    if (res.found) {
//                        chosenFaces[target_ring.half_edge().index()] = target_ring.half_edge().face().index();
//                    }
//                } while (target_ring.advance());
//
//                // construct the mass sparse matrix
//
//                /*
//                 * Prime the components
//                 * 1. interior
//                 * 2. boundary values
//                 * 3. number of boundary vertices
//                 * */
//                {
//                    for (int i = 0; i < static_cast<int>(chosenPointsIds.size()); ++i) {
//                        interior[chosenPointsIds[i]] = i;
//                        interior_rev[i] = chosenPointsIds[i];
//                    }
//
//                    for (const auto &pair: v_hat_map) {
//                        new_positions[pair.first] = pair.second;
//                    }
//
//                    num_boundary = static_cast<int>(v_hat_map.size());
//                }
//
//                /* perform math for all matrix values*/
//                {
//                    int n = static_cast<int>(chosenPoints.size());
//                    A.resize(n, n);
//                    A.setZero();
//                    Ubar = Eigen::MatrixXd::Zero(n, 1);
//                    Vbar = Eigen::MatrixXd::Zero(n, 1);
//                    for (const auto &pair: interior) {
//                        int vid = pair.first;
//                        int i = pair.second;
//                        // compute A_i
//                        {
//                            Mesh_connectivity::Vertex_ring_iterator ring = get_ring_from_target_mesh(vid);
//
//                            do {
//                                Mesh_connectivity::Vertex_iterator v_j = ring.half_edge().origin();
//
//                                // check adjacent vertex is NOT on the boundary
//                                ArraySearchResults res = arrayContains(boundary_ids, vid);
//                                if (!res.found) {
//                                    int jid = v_j.index();
//                                    int j = interior[jid];
//
//                                    // A_{ij} is inverse normalized weight for ij
//                                    A_elem.push_back(Eigen::Triplet<double>(i, j, -lambda_ij(vid, jid)));
//                                }
//                            } while (ring.advance());
//
//                        }
//
//                        // compute_UVbar_i
//                        {
//                            Mesh_connectivity::Vertex_ring_iterator ring = get_ring_from_target_mesh(vid);
//                            double u_sum = 0.0;
//                            double v_sum = 0.0;
//                            do // *points to*
//                            {
//                                Mesh_connectivity::Vertex_iterator v_j = ring.half_edge().origin();
//
//                                // check adjacent vertex IS on the boundary
//                                if (v_j.is_boundary()) {
//                                    int jid = v_j.index();
//
//                                    // normalized weight for ij
//                                    double k_ij = lambda_ij(vid, jid);
//
//                                    // new unit circle positions for j
//                                    double u_j = new_positions[jid][0];
//                                    double v_j = new_positions[jid][1];
//
//                                    double u_ij = k_ij * u_j;
//                                    double v_ij = k_ij * v_j;
//
//                                    u_sum += u_ij;
//                                    v_sum += v_ij;
//                                }
//                            } while (ring.advance());
//
//                            Ubar(i) = u_sum;
//                            Vbar(i) = v_sum;
//                        }
//                    }
//                    A.setFromTriplets(A_elem.begin(), A_elem.end());
//                }
//
//                /* solve the linear system */
//                {
//                    // Use LU since A is not symmetric!
//                    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
//                    solver.analyzePattern(A);
//                    solver.factorize(A);
//                    Eigen::MatrixXd U = solver.solve(Ubar);
//                    Eigen::MatrixXd V = solver.solve(Vbar);
//
//                    // populate vertex position map with solved UVs
//                    for (int i = 0; i < U.rows(); ++i) {
//                        new_positions[interior_rev[i]] = Eigen::Vector3d(U(i), V(i), 0.0);
//                    }
//                }
//
//                /* collect the results of the flattened section */
//                {
//                    relaxed_v_hat_map = new_positions;
//                }
//            }
//
//            /*
//             * Perform line searches to compute the optimized positions
//             * */
//            void perform_line_with_distortion_calculations() {
//                /*
//                 * Todo:
//                 *  for the target vertex,
//                 *  1. perform a line search inside the kernel v_hat
//                 *      1.1. compute the distortion as a result of this line search
//                 *      1.2. store this distortion
//                 *  2. repeat step 1 unit convergence over the PM_region as defined in Sander et. al. 2001
//                 * */
//
//                // vertex whose position we want to optimize
//                Eigen::Vector3d v = target_vertex.xyz();
//                Eigen::Vector3d optimized_v = target_vertex.xyz();
//
//                // convergence tolerance for the line search
//                double tolerance = 1.0e-5;
//
//                // original distortion
//                double o_distortion = 0.0;
//                double min_distortion = 0.0;
//
//                // distortions computed per line search
//                std::vector<double> distortions;
//
//                //
//                // begin optimization
//                //
//
//                // create list of v_hat_kernel
//                std::vector<Eigen::Vector3d> v_hat_kernel;
//                for (auto &pair: v_hat_map) {
//                    v_hat_kernel.push_back(pair.second);
//                }
//
//                double c_error = 1.0;
//                int n_iterations = 200;
//                while (c_error > tolerance) {
//                    // perform line search and compute distortion
//
//                    /*
//                     * Perform line search optimization
//                     * */
//                    {
//                        // chose a random direction of v inside v_hat_kernel
//                        int idx_1 = getRandomValue(static_cast<int>(v_hat_kernel.size()));
//                        int idx_2 = getRandomValue(static_cast<int>(v_hat_kernel.size()));
//                        Eigen::Vector3d dir = v_hat_kernel[idx_1] - v_hat_kernel[idx_2];
//                        double t = random_zero_to_one();
//                        optimized_v = optimized_v + (dir * t);
//                        double pm_error = compute_pm_error_over_target(optimized_v);
//                        c_error = pm_error - c_error;
//                    }
//
//                    /*
//                     * Compute distortion metric
//                     * */
//                    {
//                        double stretch_distortion = 0.0;
//                        double area_m1 = 0.0;
//                        double area_m2 = 0.0;
//                        double L = 0.0;
//                        double G = 0.0;
//
//                        // compute total area M1
//                        // compute total area M2
//                        // compute L and G using positions from target -> source mesh
//                        // for every triangle, compute the left and right distortion metric, and combine
//                        //   get corresponding vertex positions from map
//
//
//                        for (auto &pair: chosenFaces) {
//                            Mesh_connectivity::Face_iterator face = get_face_from_target_mesh(pair.second);
//                            Mesh_connectivity::Vertex_iterator v1 = face.half_edge().origin();
//                            Mesh_connectivity::Vertex_iterator v2 = face.half_edge().dest();
//                            Mesh_connectivity::Vertex_iterator v3 = face.half_edge().next().dest();
//                            Eigen::Vector3d q1 = v1.xyz();
//                            Eigen::Vector3d q2 = v2.xyz();
//                            Eigen::Vector3d q3 = v3.xyz();
//                            Eigen::Vector3d p1 = v_bar_map[v1.index()];
//                            Eigen::Vector3d p2 = v_bar_map[v2.index()];
//                            Eigen::Vector3d p3 = v_bar_map[v3.index()];
//
//                            double area = (((p2(0) - p1(0)) * (p3(1) - p1(1))) - ((p3(0) - p1(0)) * (p2(1) - p1(1)))) / 2;
//
//                            Eigen::Vector3d S_s = (q1 * (p2.y() - p3.y()) + q2 * (p3.y() - p1.y()) + q3 * (p1.y() - p2.y())) / 2 * area;
//                            Eigen::Vector3d S_t = (q1 * (p3.x() - p2.x()) + q2 * (p1.x() - p3.x()) + q3 * (p2.x() - p1.x())) / 2 * area;
//
//                            double a = S_s.dot(S_s);
//                            double b = S_s.dot(S_t);
//                            double c = S_t.dot(S_t);
//
//                            double L = sqrt(0.5 * ((a+c) + sqrt(pow((a-c), 2) + 4 * pow(b, 2))));
//                            double G = sqrt(0.5 * ((a+c) + sqrt(pow((a-c), 2) + 4 * pow(b, 2))));
//                        }
//
//                        // get the smallest distortion
//                        if(stretch_distortion < min_distortion) min_distortion = stretch_distortion;
//                    }
//                    n_iterations--;
//                    if (n_iterations == 0) break;
//                }
//
//                if (min_distortion >= o_distortion) {
//                    _res.state = false;
//                    _res.pos = v;
//                } else {
//                    _res.state = true;
//                    _res.pos = optimized_v;
//                }
//            }
//
//            /*
//             * Optimize a vertex.
//             * Computes the optimal position and return th
//             * */
//            OutputData optimizeVertex(const InputData &data) {
//                _opt_data = data;
//                Eigen::Vector3d optimized_pos;
//                bool isOptimized = false;
//
//                /*
//                 * todo:
//                 *      line searches with/and delaunay
//                 *      check distortion, update the results accordingly
//                 * */
//
//                // prime the optimizer
//                prime_optimizer();
//
//                // relax points to 2d space
//                relax_neighbourhood();
//
//                // line search with distortion
//                perform_line_with_distortion_calculations();
//
//                return _res;
//            }
//
//            /*
//             * Helper function
//             * ----------------
//             * Checks if a point is inside a list of points
//             * */
//
//            double compute_pm_error_over_target(Eigen::Vector3d &v) {
//                if (_opt_data.target_mesh == "m1") {
//                    return compute_pm_error_over_target(v, ism->m1());
//                } else {
//                    return compute_pm_error_over_target(v, ism->m2());
//                }
//            };
//
//            double compute_pm_error_over_target(Eigen::Vector3d &v, Mesh_connectivity &mesh) {
//                double error = 0.0;
//                for (auto &pair: v_hat_map) {
//                    error = error + error_over_neighbourhood(pair.first, pair.second, mesh);
//                }
//                return error;
//            }
//
//            double error_over_neighbourhood(int w_index, Eigen::Vector3d &v, Mesh_connectivity &mesh) {
//                double lambda = 0.5;
//                double psi = 1 / static_cast<double>(relaxed_v_hat_map.size());
//                double deviation_error = 1.0;
//                double stretch_error = 0.0;
//                double sum_area = 0.0;
//                // compute area for the flattened section
//                for (auto &pair: chosenFaces) {
//                        Mesh_connectivity::Face_iterator face = get_face_from_target_mesh(pair.second);
//                        Mesh_connectivity::Vertex_iterator v1 = face.half_edge().origin();
//                        Mesh_connectivity::Vertex_iterator v2 = face.half_edge().dest();
//                        Mesh_connectivity::Vertex_iterator v3 = face.half_edge().next().dest();
//                        Eigen::Vector3d q1 = v1.xyz();
//                        Eigen::Vector3d q2 = v2.xyz();
//                        Eigen::Vector3d q3 = v3.xyz();
//                        Eigen::Vector3d p1 = relaxed_v_hat_map[v1.index()];
//                        Eigen::Vector3d p2 = relaxed_v_hat_map[v2.index()];
//                        Eigen::Vector3d p3 = relaxed_v_hat_map[v3.index()];
//
//                        double area = (((p2(0) - p1(0)) * (p3(1) - p1(1))) - ((p3(0) - p1(0)) * (p2(1) - p1(1)))) / 2;
//
//                        Eigen::Vector3d S_s = (q1 * (p2.y() - p3.y()) + q2 * (p3.y() - p1.y()) + q3 * (p1.y() - p2.y())) / 2 * area;
//                        Eigen::Vector3d S_t = (q1 * (p3.x() - p2.x()) + q2 * (p1.x() - p3.x()) + q3 * (p2.x() - p1.x())) / 2 * area;
//
//                        double a = S_s.dot(S_s);
//                        double b = S_s.dot(S_t);
//                        double c = S_t.dot(S_t);
//
//                        double se = sqrt(0.5 * (a + c));
//
//                        sum_area = sum_area + area;
//                        stretch_error = stretch_error + se;
//                    }
//
//                return psi * (lambda * (stretch_error) + (1 - lambda) * (deviation_error * deviation_error) / sum_area);
//            }
//
//
//            bool isPointInsideBoundary(const std::vector<Eigen::Vector3d> &boundary, const Eigen::Vector3d &point) {
//                int intersectionCount = 0;
//
//                // Iterate over each edge of the boundary
//                for (size_t i = 0; i < boundary.size(); ++i) {
//                    Eigen::Vector3d p1 = boundary[i];
//                    Eigen::Vector3d p2 = boundary[(i + 1) % boundary.size()];
//
//                    // Check if the line segment intersects with the horizontal line from the point
//                    if ((p1.y() < point.y() && p2.y() >= point.y()) || (p2.y() < point.y() && p1.y() >= point.y())) {
//                        // Calculate the x-coordinate of the intersection point
//                        double intersectionX = (p1.x() + (point.y() - p1.y()) / (p2.y() - p1.y()) * (p2.x() - p1.x()));
//
//                        // Check if the intersection point is to the right of the point
//                        if (point.x() < intersectionX) {
//                            intersectionCount++;
//                        }
//                    }
//
//                    // Check if the point lies on a vertex of the boundary
//                    if ((p1.y() == point.y() && p1.x() == point.x()) || (p2.y() == point.y() && p2.x() == point.x())) {
//                        return false; // Point lies on the boundary, consider it outside
//                    }
//                }
//
//                // If the number of intersections is odd, the point is inside the boundary
//                return (intersectionCount % 2 == 1);
//            }
//
//
//            /*
//             * ==============================================
//             *              HELPER FUNCTIONS
//             * ==============================================
//             * */
//
//            /*
//             * get a random value between zero and one
//             * */
//            double random_zero_to_one() {
//                // Seed the random number generator
//                srand(time(nullptr));
//
//                // Generate a random value between 0 and RAND_MAX, and normalize it to [0, 1]
//                return static_cast<double>(rand()) / RAND_MAX;
//            }
//
//            /*
//             * returns a random interger value
//             * */
//            int getRandomValue(int n) {
//                // Seed the random number generator
//                std::random_device rd;
//                std::mt19937 gen(rd());
//
//                // Define the distribution
//                std::uniform_int_distribution<> dis(1, n);
//
//                // Generate and return a random value
//                return dis(gen);
//            }
//
//            static ArraySearchResults arrayContains(std::vector<int> &vec, int &value) {
//
//                auto it = std::find(vec.begin(), vec.end(), value);
//
//                if (it != vec.end()) {
//                    return ArraySearchResults{true, static_cast<int>(std::distance(vec.begin(), it))};
//                } else {
//                    return ArraySearchResults{false, -1};
//                }
//            }
//
//            void get_vertices_from_target(std::vector<int> &ids, std::vector<Eigen::Vector3d> &pos, int f_id, int &v1,
//                                          int &v2, int &v3) {
//
//                if (_opt_data.target_mesh == "m1") {
//
//                    Mesh_connectivity::Face_iterator f = ism->m1().face_at(f_id);
//                    Mesh_connectivity::Half_edge_iterator he = f.half_edge();
//
//                    v1 = he.origin().index();
//                    Mesh_connectivity::Vertex_iterator vv1 = ism->m1().vertex_at(v1);
//                    v2 = he.dest().index();
//                    Mesh_connectivity::Vertex_iterator vv2 = ism->m1().vertex_at(v2);
//                    v3 = he.next().dest().index();
//                    Mesh_connectivity::Vertex_iterator vv3 = ism->m1().vertex_at(v3);
//
//                    ids.push_back(v1);
//                    ids.push_back(v2);
//                    ids.push_back(v3);
//
//                    pos.push_back(vv1.xyz());
//                    pos.push_back(vv2.xyz());
//                    pos.push_back(vv3.xyz());
//
//                } else {
//
//                    Mesh_connectivity::Face_iterator f = ism->m2().face_at(f_id);
//                    Mesh_connectivity::Half_edge_iterator he = f.half_edge();
//
//                    v1 = he.origin().index();
//                    Mesh_connectivity::Vertex_iterator vv1 = ism->m2().vertex_at(v1);
//                    v2 = he.dest().index();
//                    Mesh_connectivity::Vertex_iterator vv2 = ism->m2().vertex_at(v2);
//                    v3 = he.next().dest().index();
//                    Mesh_connectivity::Vertex_iterator vv3 = ism->m2().vertex_at(v3);
//
//                    ids.push_back(v1);
//                    ids.push_back(v2);
//                    ids.push_back(v3);
//
//                    pos.push_back(vv1.xyz());
//                    pos.push_back(vv2.xyz());
//                    pos.push_back(vv3.xyz());
//                }
//            }
//
//            /*
//             * Gets the ring at the target mesh
//             * */
//            Mesh_connectivity::Vertex_ring_iterator get_ring_from_target_mesh(int v_id) {
//                if (_opt_data.target_mesh == "m1") {
//                    return ism->m1().vertex_ring_at(v_id);
//                } else {
//                    return ism->m2().vertex_ring_at(v_id);
//                }
//            }
//
//            Mesh_connectivity::Face_iterator get_face_from_target_mesh(int f_id) {
//                if (_opt_data.target_mesh == "m1") {
//                    return ism->m1().face_at(f_id);
//                } else {
//                    return ism->m2().face_at(f_id);
//                }
//            }
//
//            /*
//             * Compute lambda_ij from target
//             * */
//            double lambda_ij(int i, int j) {
//                if (_opt_data.target_mesh == "m1") {
//                    return lambda_ij_target(i, j, ism->m1());
//                } else {
//                    return lambda_ij_target(i, j, ism->m2());
//                }
//            }
//
//            /*
//             * Get the smallest value inside a vector
//             * */
//            static double get_smallest_value_in_vector(std::vector<double> &vec) {
//                // Using std::min_element
//                auto min_it = std::min_element(vec.begin(), vec.end());
//
//                if (min_it != vec.end()) {
//                    double smallest_value = *min_it;
//                    return smallest_value;
//                } else {
//                    return -1;
//                }
//
//            }
//
//            double get_wik_target(int he_index, Mesh_connectivity &mesh) {
//                // points to i from j
//                Mesh_connectivity::Half_edge_iterator he = mesh.half_edge_at(he_index);
//
//                Mesh_connectivity::Vertex_iterator I = he.dest();
//                Mesh_connectivity::Vertex_iterator J = he.origin();
//                Mesh_connectivity::Vertex_iterator P = he.twin().next().dest();
//                Mesh_connectivity::Vertex_iterator Q = he.next().dest();
//
//                Eigen::Vector3d I_xyz = I.xyz();
//                Eigen::Vector3d J_xyz = J.xyz();
//                Eigen::Vector3d P_xyz = P.xyz();
//                Eigen::Vector3d Q_xyz = Q.xyz();
//
//                Eigen::Vector3d IJ = J_xyz - I_xyz;
//                Eigen::Vector3d IP = P_xyz - I_xyz;
//                Eigen::Vector3d IQ = Q_xyz - I_xyz;
//
//                double alpha = get_angle(IJ, IP);
//                double beta = get_angle(IJ, IQ);
//                double r = IJ.norm();
//
//                return (std::tan(alpha / 2) + std::tan(beta / 2)) / r;
//
//            }
//
//            double lambda_ij_target(int i, int j, Mesh_connectivity &mesh) {
//
//                Mesh_connectivity::Vertex_iterator v_i = mesh.vertex_at(i);
//                Mesh_connectivity::Vertex_ring_iterator ring = mesh.vertex_ring_at(i);
//                double w_sum = 0.0;
//                double w_ij = 0.0;
//
//                // compute MVW for all adjacent vertices
//                do // *points to*
//                {
//                    int k = ring.half_edge().origin().index();
//                    double w_ik = get_wik_target(ring.half_edge().index(), mesh);
//
//                    // k is target adajcent vertex j
//                    if (k == j) {
//                        w_ij = w_ik;
//                    }
//
//                    w_sum += w_ik;
//
//                } while (ring.advance());
//
//                return w_ij / w_sum;
//            }
//
//            /*
//             * Computes the anngle between two vectors
//             * */
//            static double get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) {
//                return std::acos(v1.dot(v2) / (v1.norm() * v2.norm()));
//            }
//
//        private:
//
//            // the inter surface mapping this vertex belongs to
//            ISM *ism{};
//
//            // optimizer input data
//            InputData _opt_data{};
//
//            // results of the optimization process
//            OutputData _res{};
//
//            // v_bar index position map
//            // for the original target ring
//            std::map<int, Eigen::Vector3d> v_bar_map;
//
//            // v_hat index position map
//            // for flattened target ring
//            std::map<int, Eigen::Vector3d> v_hat_map;
//
//            // relaxed v_hat, v_bar positions
//            std::map<int, Eigen::Vector3d> relaxed_v_hat_map;
//
//            // v index position map
//            // for the source ring
//            std::map<int, Eigen::Vector3d> v_map;
//
//            /*
//             *
//             * Computation helpers
//             *
//             * */
//
//            Mesh_connectivity::Vertex_ring_iterator source_ring;
//            Mesh_connectivity::Vertex_ring_iterator target_ring;
//            Mesh_connectivity::Vertex_iterator target_vertex;
//            Mesh_connectivity::Vertex_iterator source_vertex;
//            std::vector<int> chosenPointsIds = {}; // v_id
//            std::vector<Eigen::Vector3d> chosenPoints = {}; // positions for chosenPoints;
//            std::map<int, int> chosenFaces = {}; // he -> face
//        };
//
//    } // end mohe
//} // end minimesh