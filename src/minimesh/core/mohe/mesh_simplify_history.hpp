
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
        };

        struct HistoryEntry {
            HistoryVertex kept_vertex_data;
            HistoryVertex removed_vertex_data;
            std::vector<HistoryHalfEdge> removed_edges;  // the first half_edge is the one between the kept and removed v
            std::vector<int> removed_faces;     // the faces belonging to the parts to be removed
        };



        class Mesh_simplify_history {
        public:
            Mesh_simplify_history(Mesh_connectivity &mesh_in) :_m(mesh_in) {}

            // Get the underlying mesh
            Mesh_connectivity & mesh() { return _m; }
            const Mesh_connectivity & mesh() const { return _m; }

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

            // creates and saves a vertex_entry
            HistoryEntry create_history_entry(int kept_v_id, int removed_v_id, std::vector<int> removed_half_edges,
                                              std::vector<int> removed_faces) {
                HistoryVertex kept_v_h = create_vertex_history(kept_v_id);
                HistoryVertex removed_v_h = create_vertex_history(removed_v_id);
                std::vector<HistoryHalfEdge> h_arr;
                h_arr.reserve(removed_half_edges.size());
                for (auto &h_id: removed_half_edges) {
                    h_arr.push_back(create_halfEdge_history(h_id));
                }
                HistoryEntry entry{kept_v_h, removed_v_h, h_arr, removed_faces};
                return entry;
            }

            // overload the create history function
            HistoryEntry create_history_entry(HistoryVertex kept_v_h, HistoryVertex removed_v_h,
                                              std::vector<HistoryHalfEdge> h_arr, std::vector<int> removed_faces) {
                HistoryEntry entry{kept_v_h, removed_v_h, h_arr, removed_faces};
                return entry;
            }

            // create and store a history
            void collect_history_entry(int kept_v_id, int removed_v_id, std::vector<int> removed_half_edges,
                                       std::vector<int> removed_faces) {
                add_to_history(create_history_entry(kept_v_id, removed_v_id, removed_half_edges, removed_faces));
            }

            // overload create and store history function
            void collect_history_entry(HistoryVertex kept_v_h, HistoryVertex removed_v_h,
                                       std::vector<HistoryHalfEdge> h_arr, std::vector<int> removed_faces) {
                add_to_history(create_history_entry(kept_v_h, removed_v_h, h_arr, removed_faces));
            }

        private:
            Mesh_connectivity & _m;
            // pointer to the mesh that we are working on.
            std::stack<HistoryEntry> history_entries;

        };
    }
}
