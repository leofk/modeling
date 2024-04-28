
#include <minimesh/core/mohe/mesh_connectivity.hpp>

namespace minimesh {

    namespace mohe {

        // structure for a history half-edge
        struct HistoryHalfEdge {
            int to;     // v_id of the origin vertex
            int from;       // v_id of the destination vertex
            int he_id;      // id of this half_edge in the original mesh. with direction `to` -> `from`
            int face_id;    // the face this half_edge belongs to
        };

        // structure for a history vertex
        struct HistoryVertex {
            int v_id;   // its index in the original mesh
            Eigen::Vector3d original_position;  // position of the vertex before contraction
            std::vector<HistoryHalfEdge> ring_edges;   // the half-edges connected to this vertex

            std::vector<int> ring_vertices() {
                std::vector<int> arr;
                for (auto &edge: ring_edges) {
                    arr.push_back(edge.from);
                }
                return arr;
            };


            // gets ring edge id
            std::vector<int> ring_edge_ids() {
                std::vector<int> arr;
                for (HistoryHalfEdge edge: ring_edges) {
                    arr.push_back(edge.he_id);
                }
                return arr;
            }
        };

        struct HistoryEntry {
            HistoryVertex kept_vertex_data;
            HistoryVertex removed_vertex_data;
            std::vector<HistoryHalfEdge> removed_edges;  // the first half_edge is the one between the kept and removed v
            std::vector<int> removed_faces;     // the faces belonging to the parts to be removed
            int top;
            int bottom;
            std::vector<int> kept_faces;
            std::vector<HistoryHalfEdge> kept_edges;

            // get kept_edges id
            std::vector<int> kept_edges_ids() {
                std::vector<int> arr;
                for (HistoryHalfEdge edge: kept_edges) {
                    arr.push_back(edge.he_id);
                }
                return arr;
            }


            // gets remove_edges id
            std::vector<int> removed_edges_ids() {
                std::vector<int> arr;
                for (HistoryHalfEdge edge: removed_edges) {
                    arr.push_back(edge.he_id);
                }
                return arr;
            }


            // gets the top and bottom points from a collapsed edge
            std::vector<int> get_adj_points() {
                std::vector<int> overlap;

                std::vector<int> ring_v1 = kept_vertex_data.ring_vertices();
                std::vector<int> ring_v2 = removed_vertex_data.ring_vertices();

                std::unordered_set<int> hash_set(ring_v1.begin(), ring_v1.end());

                // Iterate through arr2 and check for common elements
                for (int value: ring_v2) {
                    if (hash_set.count(value) > 0) {
                        overlap.push_back(value);
                    }
                }

                return overlap;
            }

        };

        /*
         * Simplification history manager,
         * Inherited by the Mesh_simplify class
         *
         * Holds edge collapse history collected at every successful edge collapse step.
         * stores this data in a stack of history entries
         *
         * A history entry, has the following
         * 1. the vertex info that was removed
         * 2. the vertex info that was kept
         * 3. top and bottom vertices adjacent to the collapsed edge
         * 4. list of deleted edges
         * 5. list of deleted faces
         * 6. list of kept faces
         * 7. list of kept edges
         * */

        class Mesh_simplify_history {
        public:
            Mesh_simplify_history(Mesh_connectivity &mesh_in) : _m(mesh_in) {}

            // Get the underlying mesh
            Mesh_connectivity &mesh() { return _m; }

            const Mesh_connectivity &mesh() const { return _m; }

            // adds an entry to the history entries
            void add_to_history(HistoryEntry entry) {
                history_entries.push(entry);
            }

            // removes and returns the most recent entry from the history
            HistoryEntry pop_history_entry() {
                HistoryEntry entry = history_entries.top();
                history_entries.pop();
                return entry;
            }

            // get the current history entries
            std::stack<HistoryEntry> history() {
                return history_entries;
            }

            // create vertex history data
            HistoryVertex create_vertex_history(int v_id) {
                Mesh_connectivity::Vertex_iterator v = _m.vertex_at(v_id);
                Mesh_connectivity::Vertex_ring_iterator ring = _m.vertex_ring_at(v_id);
                std::vector<HistoryHalfEdge> h_edges = {};

                if(v_id != v.index())printf("Create vertex error, value mismatch");

                do {
                    Mesh_connectivity::Half_edge_iterator he = ring.half_edge();
                    HistoryHalfEdge h_he{he.dest().index(), he.origin().index(), he.index(), he.face().index()};
                    h_edges.push_back(h_he);
                } while (ring.advance());

                HistoryVertex v_h{v.index(), v.xyz(), h_edges};
                return v_h;
            }

            // create half_edge history data
            HistoryHalfEdge create_halfEdge_history(int h_id) {
                Mesh_connectivity::Half_edge_iterator he = _m.half_edge_at(h_id);
                HistoryHalfEdge h{he.dest().index(), he.origin().index(), he.index(), he.face().index()};
                return h;
            }

            // creates and returns a history entry
            HistoryEntry create_history_entry(int kept_v_id, int removed_v_id, std::vector<int> removed_half_edges,
                                              std::vector<int> removed_faces, int top, int bottom,
                                              std::vector<int> kept_faces, std::vector<int> kept_edges) {
                HistoryVertex kept_v_h = create_vertex_history(kept_v_id);
                if(!mesh().vertex_at(kept_v_id).is_active())printf("Inactive vertex");
                HistoryVertex removed_v_h = create_vertex_history(removed_v_id);
                std::vector<HistoryHalfEdge> h_arr;
                std::vector<HistoryHalfEdge> kh_arr;
                h_arr.reserve(removed_half_edges.size());
                kh_arr.reserve(kept_edges.size());
                for (auto &h_id: removed_half_edges) {
                    h_arr.push_back(create_halfEdge_history(h_id));
                }
                for (auto &h_id: kept_edges) {
                    kh_arr.push_back(create_halfEdge_history(h_id));
                }
                HistoryEntry entry{kept_v_h, removed_v_h, h_arr, removed_faces, top, bottom, kept_faces, kh_arr};
                return entry;
            }

            // overload the create history function
            HistoryEntry create_history_entry(HistoryVertex kept_v_h, HistoryVertex removed_v_h,
                                              std::vector<HistoryHalfEdge> h_arr, std::vector<int> removed_faces,
                                              int top, int bottom, std::vector<int> kept_faces,
                                              std::vector<HistoryHalfEdge> kept_edges) {
                HistoryEntry entry{kept_v_h, removed_v_h, h_arr, removed_faces, top, bottom, kept_faces, kept_edges};
                return entry;
            }

            // create and store a history
            void collect_history_entry(int kept_v_id, int removed_v_id, std::vector<int> removed_half_edges,
                                       std::vector<int> removed_faces, int top, int bottom, std::vector<int> kept_faces,
                                       std::vector<int> kept_edges) {
                add_to_history(
                        create_history_entry(kept_v_id, removed_v_id, removed_half_edges, removed_faces, top, bottom,
                                             kept_faces, kept_edges));
            }

            // overload create and store history function
            void collect_history_entry(HistoryVertex kept_v_h, HistoryVertex removed_v_h,
                                       std::vector<HistoryHalfEdge> h_arr, std::vector<int> removed_faces, int top,
                                       int bottom, std::vector<int> kept_faces,
                                       std::vector<HistoryHalfEdge> kept_edges) {
                add_to_history(
                        create_history_entry(kept_v_h, removed_v_h, h_arr, removed_faces, top, bottom, kept_faces,
                                             kept_edges));
            }

        private:
            Mesh_connectivity &_m;
            // pointer to the mesh that we are working on.
            std::stack<HistoryEntry> history_entries;

        };
    }
}
