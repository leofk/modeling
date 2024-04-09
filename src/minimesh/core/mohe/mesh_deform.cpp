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
	p_prime = Eigen::MatrixXd::Zero(n-c, 3);
	Afn.resize(n,n-c);

	build_W_matrix();
	build_Afn_matrix();

    Eigen::SparseMatrix<double> AtA = Afn.transpose() * Afn;
	solver.compute(AtA);

}

void Mesh_deform::build_W_matrix() {

	int n = mesh().n_total_vertices();
	W = Eigen::MatrixXd::Zero(n, n);
	int x = 0;

	for (int i = 0; i < n; i++) {
		r_matrices[i] = Eigen::MatrixXd::Identity(3,3);
		Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
		int count = 0;
		p_all[i] = v.xyz();
		if (!is_constraint(i)) 
		{
			free[i] = x;
			free_rev[x] = i;
			p_prime.row(x) = p_all[i];
			x++;
		}	
		double sum = 0.0;

		do // *points to*
		{
			int j = ring.half_edge().origin().index();
			double wij = compute_wij(ring.half_edge().index());
			sum += wij;
			W(i, j) = wij;
			count++;
		} while (ring.advance());

		v.data().n_neighbours = count;
		W(i, i) = sum;
	} 
}

void Mesh_deform::build_Afn_matrix() {

	int n = mesh().n_total_vertices();
	
	std::vector<Eigen::Triplet<double>> Afn_elem;

	for (int i = 0; i < n; i++) {
		int mat_i = free[i]; 
		Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
		double sum = 0.0;
		
		do // *points to*
		{
			int j = ring.half_edge().origin().index();
			if (!is_constraint(j))
			{
				int mat_j = free[j];
				Afn_elem.push_back(Eigen::Triplet<double>(i, mat_j, -W(i, j)));
			} 
			sum += W(i, j);

		} while (ring.advance());

		if (!is_constraint(i)){
			Afn_elem.push_back(Eigen::Triplet<double>(i, mat_i, sum));
		}

	}

	Afn.setFromTriplets(Afn_elem.begin(), Afn_elem.end());
}

void Mesh_deform::deform(Eigen::Vector3f pull_amount) 
{

	if (!primed) {
		init();
		primed = true;
	}
	
	for(auto &pair: _handle){
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(pair.first);
		pp_handles[pair.first] = v.xyz() + pull_amount.cast<double>();
	}

	int max = 3;
	for (int i = 1; i <= max; i++) {

		Eigen::MatrixXd b = Eigen::MatrixXd::Zero(mesh().n_total_vertices(), 3);

		build_b_matrix(b);

		Eigen::MatrixXd p_prev = p_prime;
		Eigen::MatrixXd Atb = Afn.transpose() * b;
		p_prime = solver.solve(Atb);

		if (compute_energy(p_prev) <= THRESHOLD || i == max) {
			// want to save last rotation
			break;
		} else {
			r_matrices.clear();
		}
	}

	update_positions();
}


void Mesh_deform::compute_r_i(int i) 
{
	Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
	int n_neighbours = v.get_num_neighbours();

	Eigen::MatrixXd S = Eigen::MatrixXd::Zero(3, 3);
	Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);

	// Eigen::MatrixXd Pp_MAT = Eigen::MatrixXd::Zero(3, n_neighbours);
	// Eigen::MatrixXd P_MAT = Eigen::MatrixXd::Zero(3, n_neighbours);
	// Eigen::MatrixXd D = Eigen::MatrixXd::Zero(n_neighbours, n_neighbours);
	// int c = 0;

	Eigen::Vector3d ppi = p_all[i];
	Eigen::Vector3d pi = p_all[i];

	if (!is_constraint(i)) {
		ppi = p_prime.row(free[i]); 
	}

	if (_handle.count(i)>0) {
		ppi = pp_handles[i];
	}

	do // *points to*
	{
		Mesh_connectivity::Vertex_iterator vj = ring.half_edge().origin();
		int j = vj.index();

		// if j is fixed (nonhandle) constraint vertex
		Eigen::Vector3d ppj = vj.xyz();
		Eigen::Vector3d pj = p_all[j];

		if (!is_constraint(j)) {
			ppj = p_prime.row(free[j]);
		}
		if (_handle.count(j)>0) {
			ppj = pp_handles[j];
		}
	
		Eigen::Vector3d Pp = ppi - ppj;
		Eigen::Vector3d P = pi - pj; 
		S += W(i,j) * P * Pp.transpose();

		// Pp_MAT.col(c) = Pp; 
		// P_MAT.col(c) = P; 
		// D(c, c) = W(i, j);
		// c++;

	} while (ring.advance());	

	// Eigen::JacobiSVD<Eigen::MatrixXd> SVD(P_MAT * D * Pp_MAT.transpose(),
	// 										Eigen::ComputeThinU | Eigen::ComputeThinV);
	// Eigen::Matrix3d r = SVD.matrixV() * SVD.matrixU().transpose();

	// // up to changing the sign of the column of Ui corresponding
	// // to the smallest singular value, such that det (Ri) > 0.
	// if (r.determinant() < 0) {
	// 	Eigen::MatrixXd U = SVD.matrixU();
	// 	U.rightCols(1) = U.rightCols(1) * -1;
	// 	r = SVD.matrixV() * U.transpose();
	// }

	// r_matrices[i] = r;

	Eigen::JacobiSVD<Eigen::MatrixXd> SVD(S, Eigen::ComputeThinU | Eigen::ComputeThinV);

	Eigen::Matrix3d r = SVD.matrixV() * SVD.matrixU().transpose();
	Eigen::MatrixXd id = Eigen::MatrixXd::Identity(3, 3);
	id(2,2) = r.determinant();
	r_matrices[i] =  SVD.matrixV() * id * SVD.matrixU().transpose();
}

void Mesh_deform::build_b_matrix(Eigen::MatrixXd &b) {

	int n = mesh().n_total_vertices();
	
	for (int i = 0; i < n; i++) {
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);

		if (r_matrices.count(i)==0) compute_r_i(i);

		Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
		Eigen::Vector3d out = Eigen::Vector3d(0.0, 0.0, 0.0);
		
		Eigen::Vector3d pi = p_all[i]; 

		do // *points to*
		{
			Mesh_connectivity::Vertex_iterator vj = ring.half_edge().origin();
			int j = vj.index();
			if (r_matrices.count(j)==0) compute_r_i(j);

			Eigen::Vector3d pj = p_all[j];

			Eigen::Vector3d pos = pi - pj;

			Eigen::MatrixXd mul = r_matrices[i] + r_matrices[j];

			Eigen::Vector3d res = mul * pos;
			out += 0.5 * W(i,j) * res;

			if (is_constraint(j)) {
				Eigen::Vector3d ppj = p_all[j];
				if (_handle.count(j)>0) {
					ppj = pp_handles[j];
				}
				out += W(i,j) * ppj;
			}

		} while (ring.advance());

		if (is_constraint(i)) {
			Eigen::Vector3d ppi = p_all[i];
			if (_handle.count(i)>0) {
				ppi = pp_handles[i];
			}
			out -= W(i,i) * ppi;
		}

		b.row(i) = out;

	}
}

bool Mesh_deform::is_constraint(int vid) {
	return _constraints.count(vid) > 0;
}

double Mesh_deform::compute_energy(Eigen::MatrixXd &p_prev) {

	// Initialize sum of squared distances
    double totalSquaredDistance = 0.0;

    // Compute the sum of squared distances between corresponding points
    for (int i = 0; i < p_prime.rows(); ++i) {
        // Compute the squared Euclidean distance between corresponding points
        double squaredDistance = (p_prime.row(i) - p_prev.row(i)).squaredNorm();

        // Add the squared distance to the total
        totalSquaredDistance += squaredDistance;
    }

    // Compute the mean squared distance
    double meanSquaredDistance = totalSquaredDistance / p_prime.rows();

    // Compute the mean distance by taking the square root of the mean squared distance
    double meanDistance = sqrt(meanSquaredDistance);

    return meanDistance;
}

void Mesh_deform::update_positions() {
	for (int i = 0; i < mesh().n_total_vertices(); ++i) {
		if (!is_constraint(i)) 
		{
			Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
			v.data().xyz = p_prime.row(free[i]);
		}
	}

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

// is boundary edge
bool Mesh_deform::is_boundary_edge(int he_index) {
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_index);

	// see if the face for given half edge (or its twin) has negative index
	return he.twin().face().index() < 0 || he.face().index() < 0;
}

//
// compute cotangent weights for a half edge
//
double Mesh_deform::compute_wij(int he_id) {
	// points to i from j
	Mesh_connectivity::Half_edge_iterator he = mesh().half_edge_at(he_id);
            // if we are on the boundary
	if (is_boundary_edge(he_id)) {

		// get the edge inside the mesh
		if (he.face().index() < 0) he = he.twin();

		Mesh_connectivity::Vertex_iterator p1 = he.dest();
		Mesh_connectivity::Vertex_iterator p2 = he.origin();
		Mesh_connectivity::Vertex_iterator p3 = he.next().dest();
		
		Eigen::Vector3d p3_xyz = p3.xyz();
		// if (!is_constraint(p3.index())) {
		// 	p3_xyz = p_free[p3.index()]; 
		// }

		Eigen::Vector3d p2_xyz = p2.xyz();
		// if (!is_constraint(p2.index())) {
		// 	p2_xyz = p_free[p2.index()]; 
		// }
		Eigen::Vector3d p1_xyz = p1.xyz();
		// if (!is_constraint(p1.index())) {
		// 	p1_xyz = p_free[p1.index()]; 
		// }

		Eigen::Vector3d v1 = p3_xyz - p1_xyz;
		Eigen::Vector3d v2 = p3_xyz - p2_xyz;

		double alpha = get_angle(v1, v2);
		double cot_alpha = 0.5 * (1 / std::tan(alpha));
		return cot_alpha;
	}

	Mesh_connectivity::Vertex_iterator I = he.dest();
	Mesh_connectivity::Vertex_iterator J = he.origin();
	Mesh_connectivity::Vertex_iterator P = he.twin().next().dest();
	Mesh_connectivity::Vertex_iterator Q = he.next().dest();

	Eigen::Vector3d I_xyz = I.xyz();
	// if (!is_constraint(I.index())) {
	// 	I_xyz = p_free[I.index()]; 
	// }
	Eigen::Vector3d J_xyz = J.xyz();
	// if (!is_constraint(J.index())) {
	// 	J_xyz = p_free[J.index()]; 
	// }
	Eigen::Vector3d P_xyz = P.xyz();
	// if (!is_constraint(P.index())) {
	// 	P_xyz = p_free[P.index()]; 
	// }
	Eigen::Vector3d Q_xyz = Q.xyz();
	// if (!is_constraint(Q.index())) {
	// 	Q_xyz = p_free[Q.index()]; 
	// }

	Eigen::Vector3d IP = P_xyz - I_xyz;
	Eigen::Vector3d JP = P_xyz - J_xyz;
	Eigen::Vector3d IQ = Q_xyz - I_xyz;
	Eigen::Vector3d JQ = Q_xyz - J_xyz;

	double alpha = get_angle(IQ, JQ);
	double beta = get_angle(IP, JP);

	return 0.5 * ((1 / std::tan(alpha)) + (1 / std::tan(beta)));
	// return 1.0;
}

} // end of mohe
} // end of minimesh/