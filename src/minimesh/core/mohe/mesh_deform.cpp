#include <minimesh/core/mohe/mesh_deform.hpp>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <omp.h>

namespace minimesh {
namespace mohe {

void Mesh_deform::init() {
	int n = mesh().n_total_vertices();
	int c = _constraints.size();
	bf = Eigen::MatrixXd::Zero(n-c, 3);
	xc = Eigen::MatrixXd::Zero(c, 3);
	Afc = Eigen::MatrixXd::Zero(n-c, c);

	// Aff = Eigen::MatrixXd::Zero(n-c, n-c);
    Eigen::SparseMatrix<double> Aff;
	Aff.resize(n-c,n-c);

	build_c_map();
	build_W_matrix();
	build_L_matrix(Aff);

	solver.compute(Aff);

}

void Mesh_deform::build_c_map() 
{
	int i = 0;
	for (const auto& pair : _constraints) 
	{
		c_map[pair.first] = i;
		i++;
	}
}

void Mesh_deform::build_W_matrix() {
	// build nxn mat
	// 1 ; when i==j
	// wij ; when j in N(i)
	// 0 otherwise

	int n = mesh().n_total_vertices();
	W = Eigen::MatrixXd::Zero(n, n);

	// p_prev = Eigen::MatrixXd::Zero(n, 3);
	int x = 0;

	for (int i = 0; i < n; i++) {
		Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
		int count = 0;

		if (!is_constraint(i)) 
		{
			free[i] = x;
			free_rev[x] = i;
			x++;
		}	

		do // *points to*
		{
			int j = ring.half_edge().origin().index();
			W(i, j) = compute_wij(ring.half_edge().index());
			count++;
		} while (ring.advance());

		v.data().n_neighbours = count;
		W(i, i) = 1.0;

		// I_matrices[i] = Eigen::MatrixXd::Identity(3, 3);
		r_matrices[i] = Eigen::MatrixXd::Identity(3, 3);
		// p_prev.row(i) = v.xyz();
	} 
}

void Mesh_deform::build_L_matrix(Eigen::SparseMatrix<double> &Aff) {
	// build nxn mat
	// sum wij for j in N(i) ; when i==j
	// -wij ; when j in N(i)
	// 0 ; otherwise

	int n = mesh().n_total_vertices();
	
	std::vector<Eigen::Triplet<double>> Aff_elem;


	for (int i = 0; i < n; i++) {
		if (!is_constraint(i))
		{
			int mat_i = free[i]; 
			Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
			double sum = 0.0;
			
			do // *points to*
			{
				int j = ring.half_edge().origin().index();
				if (!is_constraint(j))
				{
					int mat_j = free[j];
					sum += W(i, j);
					Aff_elem.push_back(Eigen::Triplet<double>(mat_i, mat_j, -W(i, j)));
					// Aff(mat_i,mat_j) = -W(i, j);
				} else {
					int mat_j = c_map[j];
					sum += W(i, j);
					Afc(mat_i,mat_j) = -W(i, j);
				}

			} while (ring.advance());

			Aff_elem.push_back(Eigen::Triplet<double>(mat_i, mat_i, sum));
			// Aff(mat_i,mat_i) = sum;
		} else {
		// if the item is a handle.
		if(_handle.count(i)>0){
			pp_handles[i] = mesh().vertex_at(i).xyz();
		}
	}
	}

	Aff.setFromTriplets(Aff_elem.begin(), Aff_elem.end());
}

void Mesh_deform::deform(int _handle_id, Eigen::Vector3f pull_amount) 
{

	if (!primed) {
		init();

		primed = true;
		printf("primed");
	}

	handle_id = _handle_id;
	Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(handle_id);
	// pp_handle = v.xyz() + pull_amount.cast<double>();
	// shift all the handles by the pull amount
	for(auto &pair: pp_handles){
		pp_handles[pair.first] = pair.second + pull_amount.cast<double>();
	}

	// while not converged
	// compute rotation mats per vertex (when i>0)
	// compute b matrix
	// solve for p'
	// update positions
	// compute energy
	bool converged = false;
	bool first = true;

	for (int i = 0; i < 1; i++) {
		if (!first) {
			// use identify matrices in R for first iteration
			compute_r_matrices();
			// printf("r \n");
		} else {
			// r_matrices = I_matrices;
			first = false;
		}

		build_b_matrix();

		Eigen::MatrixXd b = bf - Afc * xc;
		p_prime = solver.solve(b);
	}
	xc_done = false;

	update_positions();
}


void Mesh_deform::compute_r_matrices() {
	// for each vertex
	// setup P, P_prime, D
	// compute SVD
	// r = VU^T
	int n = mesh().n_total_vertices();

	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
		int n_neighbours = v.get_num_neighbours();

		Eigen::MatrixXd Pp = Eigen::MatrixXd::Zero(3, n_neighbours);
		Eigen::MatrixXd P = Eigen::MatrixXd::Zero(3, n_neighbours);
		Eigen::MatrixXd D = Eigen::MatrixXd::Zero(n_neighbours, n_neighbours);

		Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
		
		Eigen::Vector3d ppi = v.xyz();
		if (!is_constraint(i)) ppi = p_prime.row(free[i]);
		// if (i==handle_id) ppi = pp_handle;
		if (pp_handles.count(i)>0) ppi = pp_handles[i];

		int c = 0;
		do // *points to*
		{
			Mesh_connectivity::Vertex_iterator vj = ring.half_edge().origin();
			int j = vj.index();

			Eigen::Vector3d ppj = vj.xyz();
			if (!is_constraint(j)) ppj = p_prime.row(free[j]);
			// if (j==handle_id) ppj = pp_handle;
			if (pp_handles.count(j)>0) ppj = pp_handles[j];

			Pp.col(c) = ppi - ppj;; // this should be new pos
			P.col(c) = v.xyz() - vj.xyz(); // this should be old pos
			D(c, c) = W(i, j);
			c++;
		} while (ring.advance());

		Eigen::JacobiSVD<Eigen::MatrixXd> SVD(P * D * Pp.transpose(),
												Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::Matrix3d r = SVD.matrixV() * SVD.matrixU().transpose();

		// up to changing the sign of the column of Ui corresponding
		// to the smallest singular value, such that det (Ri) > 0.
		if (r.determinant() < 0) {
			Eigen::MatrixXd U = SVD.matrixU();
			U.rightCols(1) = U.rightCols(1) * -1;
			r = SVD.matrixV() * U.transpose();
		}

		r_matrices[i] = r;
	}
}

void Mesh_deform::build_b_matrix() {
	// for each vertex
	// if v is a constraint
	// if fixed = pos is orig
	// if handle = pos is new
	// else
	// compute RHS of (8)

	int n = mesh().n_total_vertices();

	for (int i = 0; i < n; i++) {
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);

		if (is_constraint(i) && !xc_done) {
			// TODO IS THIS THE CASE FOR THE HANDLE? OR SHOULD WE HANDLE THE HANDLE INDIV
			// if not, make sure we update the position when handle if first assignment in modifier
			Eigen::Vector3d res = v.xyz();
			// if (i==handle_id) res = pp_handle;
			if (pp_handles.count(i)>0) res = pp_handles[i];
			xc.row(c_map[i]) = res;
		} 

		if (!is_constraint(i))
		{
			Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
			Eigen::Vector3d out = Eigen::Vector3d(0.0, 0.0, 0.0);

			do // *points to*
			{
				Mesh_connectivity::Vertex_iterator vj = ring.half_edge().origin();
				int j = vj.index();
				Eigen::Vector3d pos = v.xyz() - vj.xyz();
				Eigen::MatrixXd mul = r_matrices[i] + r_matrices[j];
				Eigen::Vector3d res = mul * pos;
				out += 0.5 * W(i,j) * res;

			} while (ring.advance());
			bf.row(free[i]) = out;
		}
	}
	xc_done = true;
}

bool Mesh_deform::is_constraint(int vid) {
	// for (int c: constraint_vids) {
	// 	if (c == vid) return true;
	// }
	// return false;
	return _constraints.count(vid) > 0;
}

double Mesh_deform::compute_energy() {
	double energy = 0.0;
	int n = mesh().n_total_vertices();

	for (int i = 0; i < n; i++) {

		Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
		double cell_energy = 0.0;

		do // *points to*
		{
			Mesh_connectivity::Vertex_iterator vj = ring.half_edge().origin();
			int j = vj.index();
			// todo - below will break down for pprime when i is constraint
			Eigen::Vector3d a1 = p_prime.row(free[i]) - p_prime.row(free[j]);
			Eigen::Vector3d a2 = v.xyz() - vj.xyz();
			Eigen::Vector3d mul = r_matrices[i] * a2;
			Eigen::Vector3d v = (a1) - mul;
			cell_energy += W(i, j) * v.squaredNorm();
		} while (ring.advance());

		energy += W(i, i) * cell_energy;
	}
	return energy;
}

void Mesh_deform::update_positions() {
	// update vertex positions
	for (int i = 0; i < mesh().n_total_vertices(); ++i) {
		if (!is_constraint(i)) 
		{
			Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
			v.data().xyz = p_prime.row(free[i]);
		}
	}

	// Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(handle_id);
	// v.data().xyz = pp_handle;

	for(auto &pair: pp_handles)
	{
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(pair.first);
		v.data().xyz = pair.second;
	}
}

// 
// compute angle between two vectors (cosine similarity)
//
double Mesh_deform::get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) {
	return std::acos(v1.dot(v2) / (v1.norm() * v2.norm()));
}

//
// compute cotangent weights for a half edge
//
double Mesh_deform::compute_wij(int he_id) {
	// points to i from j
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);

	Mesh_connectivity::Vertex_iterator I = he.dest();
	Mesh_connectivity::Vertex_iterator J = he.origin();
	Mesh_connectivity::Vertex_iterator P = he.twin().next().dest();
	Mesh_connectivity::Vertex_iterator Q = he.next().dest();

	Eigen::Vector3d I_xyz = I.xyz();
	Eigen::Vector3d J_xyz = J.xyz();
	Eigen::Vector3d P_xyz = P.xyz();
	Eigen::Vector3d Q_xyz = Q.xyz();

	Eigen::Vector3d IP = P_xyz - I_xyz;
	Eigen::Vector3d JP = P_xyz - J_xyz;
	Eigen::Vector3d IQ = Q_xyz - I_xyz;
	Eigen::Vector3d JQ = Q_xyz - J_xyz;

	double alpha = get_angle(IQ, JQ);
	double beta = get_angle(IP, JP);

	return 0.5 * ((1 / std::tan(alpha)) + (1 / std::tan(beta)));
}

} // end of mohe
} // end of minimesh