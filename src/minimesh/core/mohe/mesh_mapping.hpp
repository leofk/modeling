#pragma once

#include <string>

#include <minimesh/core/mohe/mesh_connectivity.hpp>
#include <minimesh/core/mohe/mesh_simplify.hpp>
#include <minimesh/core/util/numbers.hpp>
#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <string>
#include <random>
#include <unordered_set>
#include <ctime>
#include <cmath>
#include <iostream>

#define M1 1
#define M2 2

struct Inter_data {
    std::vector<double> bary;
    int f_id;
};

namespace minimesh {
    namespace mohe {
        /*
         * Structure of data for calling the optimization function in the optimizer class
         * Source mesh is the mesh we are coming from, with default value M1
         * Target mesh it the mesh we are going to, with default value M2
         * v_id is the vertex whose position we want to optimize from m2 to m1
         *
         * */
        struct InputData {
            int source_mesh = M2;
            int target_mesh = M1;
            int source_v_id;
        };

        struct OutputData {
            bool state = true;
            Eigen::Vector3d pos = Eigen::Vector3d(0, 0, 0);
        };

        struct ArraySearchResults {
            bool found = false;
            int index = -1;
        };

        class Mesh_mapping {
        public:

            int OPT_COUNT = 0;

            /*
             * ================================================================
             * VERTEX OPTIMIZER
             * -----------------------------------------------------------------
             * Class definition.
             * Forward declaration and then skeletal definitions
             * =================================================================
             * */

            class VertexOptimizer;

            class VertexOptimizer {
            public:

                /* Get the optimization results*/
                OutputData res() { return _res; }

                /*
                 * Gets the ring for the target mesh, computes the flattened positions in 2d space and saves them their
                 * required mapped regions.
                 * */
                void prime_optimizer();

                /*
                 * Mean value parameterization to relax points in N(v_bar) into N(v_hat)
                 * */
                void relax_neighbourhood();

                /*
                 * Perform line searches to compute the optimized positions
                 * */
                void perform_line_with_distortion_calculations();


                /*
                 * Optimize a vertex.
                 * Computes the optimal position and return the results if the system was optimized
                 * */
                OutputData optimizeVertex(const InputData &data);


                /*
                 * ==============================================
                 *        HELPER FUNCTIONS FOR OPTIMIZER
                 * ==============================================
                 * */

                /*
                 * Gets the image position from the parent
                 * */
                Eigen::Vector3d get_image_position(int v_id);

                double compute_pm_error_over_target(Eigen::Vector3d &v);

                double compute_pm_error_over_target(Eigen::Vector3d &v, Mesh_connectivity &mesh);

                double error_over_neighbourhood(int w_index, Eigen::Vector3d &v, Mesh_connectivity &mesh);

                bool
                isPointInsideBoundary(const std::vector<Eigen::Vector3d> &boundary, const Eigen::Vector3d &point);

                /*
                 * get a random value between zero and one
                 * */
                static double random_zero_to_one() {
                    // Seed the random number generator
                    srand(time(nullptr));

                    // Generate a random value between 0 and RAND_MAX, and normalize it to [0, 1]
                    return static_cast<double>(rand()) / RAND_MAX;
                }

                /*
                 * returns a random interger value
                 * */
                static int getRandomValue(int n) {
                    // Seed the random number generator
                    std::random_device rd;
                    std::mt19937 gen(rd());

                    // Define the distribution
                    std::uniform_int_distribution<> dis(1, n);

                    // Generate and return a random value
                    return dis(gen);
                }

                static ArraySearchResults arrayContains(std::vector<int> &vec, int &value) {

                    auto it = std::find(vec.begin(), vec.end(), value);

                    if (it != vec.end()) {
                        return ArraySearchResults{true, static_cast<int>(std::distance(vec.begin(), it))};
                    } else {
                        return ArraySearchResults{false, -1};
                    }
                }


                /*
                 * Gets the ring at the target mesh
                 * */
                Mesh_connectivity::Vertex_ring_iterator get_ring_from_target_mesh(int v_id);

                Mesh_connectivity::Face_iterator get_face_from_target_mesh(int f_id);

                /*
                 * Compute lambda_ij from target
                 * */
                double lambda_ij(int i, int j);

                double get_wik_target(int he_index, Mesh_connectivity &mesh);

                double lambda_ij_target(int i, int j, Mesh_connectivity &mesh);

                Mesh_connectivity::Vertex_ring_iterator ring_at_source();

                /*
                 * Computes the angle between two vectors
                 * */
                static double get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) {
                    return std::acos(v1.dot(v2) / (v1.norm() * v2.norm()));
                }

                /* default constructor */
                VertexOptimizer() =  default;

            private:

                // the inter surface mapping this vertex belongs to
                Mesh_mapping * ism{};

                // optimizer input data
                InputData _opt_data;

                // results of the optimization process
                OutputData _res;

                // v_bar index position map
                // for the original target ring
                std::map<int, Eigen::Vector3d> v_bar_map;

                // v_hat index position map
                // for flattened target ring
                std::map<int, Eigen::Vector3d> v_hat_map;

                // relaxed v_hat, v_bar positions
                std::map<int, Eigen::Vector3d> relaxed_v_hat_map;

                // v index position map
                // for the source ring
                std::map<int, Eigen::Vector3d> v_map;

                /*
                 *
                 * Computation helpers
                 *
                 * */

                Mesh_connectivity::Vertex_iterator source_vertex;
                std::vector<int> chosenPointsIds; // v_id
                std::vector<Eigen::Vector3d> chosenPoints; // positions for chosenPoints;
                std::map<int, int> chosenFaces; // he -> face

                VertexOptimizer(Mesh_mapping *meshMapping): ism(meshMapping){}
                friend class Mesh_mapping;
            };



            /*
             * ================================================================
             * MESH MAPPING CONTINUES
             * =================================================================
             * */


            Mesh_mapping(Mesh_connectivity &m1, Mesh_connectivity &m2) : _m1(m1), _m2(m2) {}

            Mesh_connectivity &mesh() {
                if (current_mesh == M1)
                    return _m1;
                else
                    return _m2;
            }

            void switch_mesh() {
                if (current_mesh == M1)
                    current_mesh = M2;
                else
                    current_mesh = M1;
            }

            std::map<int, Inter_data> &current_ISM() {
                if (current_mesh == M1)
                    return ISM_M1;
                else
                    return ISM_M2;
            }

            std::map<int, int> &current_FTV() {
                if (current_mesh == M1)
                    return m1_ftv_map;
                else
                    return m2_ftv_map;
            }

            void build_mapping();

            void init_ISM();

            Inter_data ISM_iteration(HistoryEntry &entry);

            Inter_data trivial_map_data(int vid);

            Inter_data find_enclosing_face(std::vector<int> faces, Eigen::Vector3d v_pos);

            Eigen::Vector3d get_pos_from_inter_data(Inter_data data);

            std::vector<double> compute_bc(int f_id, Eigen::Vector3d v_pos);

            double get_wij(Eigen::Vector3d i, Eigen::Vector3d j, std::vector<Eigen::Vector3d> &verts);

            void update_positions(int mesh);

            void split_vertex(Mesh_simplify &simp_m);

            double get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);

            // return M2
            Mesh_connectivity &m1() { return _m1; }

            // return M1
            Mesh_connectivity &m2() { return _m2; }

            // return a vertex optimizer
            VertexOptimizer vertexOptimizer();

        private:
            Mesh_connectivity &_m1;
            Mesh_connectivity &_m2;

            std::map<int, int> mesh_map = {
                    // cow1 to horse init map
                    {2348, 6843},
                    {2149, 9676},
                    {1782, 9078},
                    {1450, 5922},
            };

            std::map<int, int> m1_ftv_map; // M1 Face to M2 Vertex Map
            std::map<int, int> m2_ftv_map; // M2 Face to M1 Vertex Map
            std::map<int, Inter_data> ISM_M1; // M1 Inter-surface Map
            std::map<int, Inter_data> ISM_M2; // M2 Inter-surface Map
            int current_mesh; // Current working mesh
        };


    } // end of mohe
} // end of minimesh
