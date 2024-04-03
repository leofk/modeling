#include <minimesh/core/mohe/Mesh_deform.hpp>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace minimesh
{
namespace mohe
{

void Mesh_deform::init()
{
	int n = mesh().n_total_vertices();
    b = Eigen::MatrixXd::Zero(n, 3);   


	build_W_matrix();
	build_L_matrix();
}


void Mesh_deform::build_W_matrix() 
{
	// build nxn mat
	// 1 ; when i==j
	// wij ; when j in N(i)
	// 0 otherwise

	int n = mesh().n_total_vertices();
	W = Eigen::MatrixXd::Zero(n, n);
    
	p_prev = Eigen::MatrixXd::Zero(n, 3);   

    for (int i = 0; i < n; i++) 
	{
		Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
		int count;

		do // *points to*
		{
			int j = ring.half_edge().origin().index();
			W(i,j) = compute_wij(ring.half_edge().index());
			count++;
		} while(ring.advance());

		v.data().n_neighbours = count;
		W(i,i) = 1.0;

        r_matrices[i] = Eigen::MatrixXd::Identity(3, 3);
        p_prev.row(i) = v.xyz();	
	}
}

void Mesh_deform::build_L_matrix() 
{
	// build nxn mat
	// sum wij for j in N(i) ; when i==j
	// -wij ; when j in N(i)
	// 0 ; otherwise

	int n = mesh().n_total_vertices();

    for (int i = 0; i < n; i++) 
	{
		Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);
		double sum;
		do // *points to*
		{
			int j = ring.half_edge().origin().index();
			sum += W[i,j];
			L_elem.push_back(Eigen::Triplet<double>(i, j, -W[i,j]));

		} while(ring.advance());

		L_elem.push_back(Eigen::Triplet<double>(i, i, sum));
	}

	L.setFromTriplets(L_elem.begin(), L_elem.end());
}

void Mesh_deform::deform(std::vector<int> _constraint_vids, int _handle_vid) 
{
	constraint_vids.clear();
	constraint_vids = _constraint_vids;
	handle_vid = _handle_vid;

	update_L_matrix();

	// while not converged
		// compute rotation mats per vertex (when i>0)
		// compute b matrix
		// solve for p'
		// update positions
		// compute energy
	bool converged = false;
	bool first = true;

	int i;
	while(!converged)
	{
		if (!first) { 
			// use identify for first iteration
			compute_r_matrices();
		} else {
			first = false;
		}


		build_b_matrix();

		p_prev = p_prime;
		
		solve_system();

		// TODO perhaps below can be just do done after i iterations?
		// if (compute_energy() <= THRESHOLD) converged = true;
		if (i==10) converged = true;
		i++;
	}

	update_positions();
}

void Mesh_deform::compute_r_matrices() 
{
	// for each vertex
		// setup P, P_prime, D
		// compute SVD
		// r = VU^T
	int n = mesh().n_total_vertices();

    for (int i = 0; i < n; i++) 
	{
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
        int n_neighbours = v.get_num_neighbours();

        Eigen::MatrixXd Pp = Eigen::MatrixXd::Zero(3, n_neighbours);
        Eigen::MatrixXd P = Eigen::MatrixXd::Zero(3, n_neighbours);
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(n_neighbours, n_neighbours);

		Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);

		do // *points to*
		{
			int j = ring.half_edge().origin().index();
            Pp.col(j) = p_prime.row(i) - p_prime.row(j); // this should be new pos
            P.col(j) = p_prev.row(i) - p_prev.row(j); // this should be old pos
            D(j, j) = W(i, j);
		} while(ring.advance());

        Eigen::JacobiSVD<Eigen::MatrixXd> SVD(P * D * Pp.transpose());
       	Eigen::Matrix3d r = SVD.matrixV() * SVD.matrixU().transpose(); 
        
		// up to changing the sign of the column of Ui corresponding 
		// to the smallest singular value, such that det (Ri) > 0.
		if(r.determinant() < 0)
        {
            Eigen::MatrixXd U = SVD.matrixU();
            U.rightCols(1) = U.rightCols(1) * -1; 
            r = SVD.matrixV() * U.transpose(); 
        }

		r_matrices[i] = r;
    }
}

void Mesh_deform::build_b_matrix() 
{
	// for each vertex
		// if v is a constraint
			// if fixed = pos is orig
			// if handle = pos is new
		// else
			// compute RHS of (8)

	int n = mesh().n_total_vertices();

	for ( int i = 0; i < n; i++)
    {
        Eigen::Vector3d pos;
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
		
		if(is_constraint(i))
		{
			// TODO IS THIS THE CASE FOR THE HANDLE? OR SHOULD WE HANDLE THE HANDLE INDIV
			// if not, make sure we update the position when handle if first assignment in modifier
			pos = v.xyz(); 
        }
        else
		{
        	Eigen::Vector3d sum;
			Mesh_connectivity::Vertex_ring_iterator ring = mesh().vertex_ring_at(i);

			do // *points to*
			{
				int j = ring.half_edge().origin().index();
				Mesh_connectivity::Vertex_iterator v_j = mesh().vertex_at(j);
                sum += W(i,j) * 0.5 * (r_matrices[i] + r_matrices[j]) * (v.xyz() - v_j.xyz());

			} while(ring.advance());

			pos = sum;
        }

        b.row(i) = pos;
    }
}

bool Mesh_deform::is_constraint(int vid) 
{
	for (int c : constraint_vids) {
		if (c == vid) return true;
	}
	return false;
}

void Mesh_deform::solve_system() 
{
	// solve L*p_prime = b using sparse cholesky solver

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    // Compute the symbolic decomposition
    solver.analyzePattern(L);

    // Compute the numeric decomposition
    solver.factorize(L);

    // Solve the system
    p_prime = solver.solve(b);
}

double Mesh_deform::compute_energy()
{
	return 0.0;
}

void Mesh_deform::update_positions()
{
	// update vertex positions
	for(int i = 0 ; i < mesh().n_total_vertices() ; ++i)
	{
		Mesh_connectivity::Vertex_iterator v = mesh().vertex_at(i);
		v.data().xyz = p_prime.row(i);
	}
}

void Mesh_deform::update_L_matrix()
{
	// for i in constraint_vids
		// Lij =
		// 1 ; when i==j
		// 0 ; otherwise	
	
    for (int i = 0; i < constraint_vids.size(); i++) 
	{
		int vid = constraint_vids[i];
        // Modify L(i, j) for each row i
        for (int j = 0; j < L.cols(); j++) 
        {
            if (vid == j) {
                // Set L(i, j) to 1 if i == j
                L.coeffRef(vid, j) = 1.0;
            } else {
                // Set L(i, j) to 0 otherwise
                L.coeffRef(vid, j) = 0.0;
            }
        }
	}
}

// 
// compute angle between two vectors (cosine similarity)
//
double Mesh_deform::get_angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) 
{
	return std::acos(v1.dot(v2) / (v1.norm() * v2.norm()));
}

//
// compute cotangent weights for a half edge
//
double Mesh_deform::compute_wij(int he_id) 
{
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

	Eigen::Vector3d IP = P_xyz-I_xyz;
	Eigen::Vector3d JP = P_xyz-J_xyz;
	Eigen::Vector3d IQ = Q_xyz-I_xyz;
	Eigen::Vector3d JQ = Q_xyz-J_xyz;

	double alpha = get_angle(IQ, JQ);
	double beta = get_angle(IP, JP);

    return 0.5 * ((1 / std::tan(alpha)) + (1 / std::tan(beta)));
}

} // end of mohe
} // end of minimesh
