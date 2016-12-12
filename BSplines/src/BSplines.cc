#include <iostream>
#include <MABSR/BSplines.h>
#include <Eigen/Sparse>

int verbose;

namespace mabsr {
	BSplines::BSplines() {}

	BSplines::BSplines(const Vector &_knots) : knots_(_knots) {}

	Scalar BSplines::B(const int &i, const int &k, const Scalar &t) const {
		if(t >= knots_(i+k) || t < knots_(i)) return 0.0;
		if(k == 1) return X(i, t);
		return w(i, k, t) * B(i, k-1, t) + (1 - w(i+1, k, t)) * B(i+1, k-1, t);
	}

	Scalar BSplines::dB(const int &i, const int &k, const Scalar &t) const {
		if(k <= 1) return 0.0;

		Scalar rt1 = 0.0, rt2 = 0.0;
		if(knots_(i+k) != knots_(i+1)) rt1 = -B(i+1, k-1, t) / (knots_(i+k) - knots_(i+1));
		if(knots_(i+k-1) != knots_(i)) rt2 = B(i, k-1, t) / (knots_(i+k-1) - knots_(i));

		return (k-1) * (rt1 + rt2);
	}

	Scalar BSplines::d2B(const int &i, const int &k, const Scalar &t) const {
		if(k <= 2) return 0.0;

		Scalar rt1 = 0.0, rt2 = 0.0;
		if(knots_(i+k) != knots_(i+1)) rt1 = -dB(i+1, k-1, t) / (knots_(i+k) - knots_(i+1));
		if(knots_(i+k-1) != knots_(i)) rt2 = dB(i, k-1, t) / (knots_(i+k-1) - knots_(i));

		return (k-1) * ( rt1 + rt2);
	}

	Scalar BSplines::w(const int &i, const int &k, const Scalar &t) const {
		if (knots_(i) != knots_(i+k-1)) {
			return (t - knots_(i)) / (knots_(i+k-1) - knots_(i));
		}
		return 0.0;
	}

	Scalar BSplines::X(const int &i, const Scalar &t) const {
		if(t >= knots_(i) && t < knots_(i+1)) return 1;
		return 0.0;
	}

	Vector equally_spaced_bar_t(const PointSet &q) {
		int M = q.rows() - 1;
		Vector bar_t(M+1);
		for(int j = 0; j < M+1; j++) bar_t(j) = (Scalar)j/M;
		return bar_t;
	}

	Vector chordal_bar_t(const PointSet &q) {
		int M = q.rows() - 1;
		Vector bar_t(M+1);
		Scalar d = 0.0;
		for(int j = 1; j < M+1; j++) d += (q.row(j) - q.row(j-1)).norm();
		bar_t(0) = 0.0;
		for(int j = 1; j < M+1; j++) bar_t(j) =	bar_t(j-1) + (q.row(j) - q.row(j-1)).norm() / d;
		bar_t(M) = 1.0;
		return bar_t;
	}

	Vector centripetal_bar_t(const PointSet &q) {
		int M = q.rows() - 1;
		Vector bar_t(M+1);
		Scalar d = 0.0;
		for(int j = 1; j < M+1; j++) d += sqrt( (q.row(j) - q.row(j-1)).norm() );
		bar_t(0) = 0.0;
		for(int j = 1; j < M+1; j++) bar_t(j) =	bar_t(j-1) + sqrt( (q.row(j) - q.row(j-1)).norm() ) / d;
		bar_t(M) = 1.0;
		return bar_t;
	}

	void adapt_bar_t(Vector &bar_t, const Scalar &scale, const Scalar &shift) {
		bar_t *= scale;
		bar_t.array() += shift;
	}

	Vector equally_spaced_knots(const int &M, const int &k, const bool &natural_end) {
		int n = k-1;
		if (natural_end) {
			Vector knots(M+3+k);
			for (int i = 0; i < n+1; i++) knots(i) = 0;
			for (int i = M+3; i < M+k+3; i++) knots(i) = 1;
			for(int i = 1; i < M-n+3; i++) knots(n+i) = (Scalar)i / (M+3-n);
			return knots;

		} else {
			Vector knots(M+1+k);
			for (int i = 0; i < n+1; i++) knots(i) = 0;
			for (int i = M+1; i < M+k+1; i++) knots(i) = 1;
			for(int i = 1; i < M-n+1; i++) knots(n+i) = (Scalar)i / (M+1-n);
			return knots;
		}
	}

	Vector average_knots(const int &M, const int &k, const Vector &bar_t, const bool &natural_end) {
		assert(bar_t.size() == (M+1));

		if (natural_end) {
			int n = k-1;
			Vector knots(M+3+k);
			for (int i = 0; i < n+1; i++) knots(i) = 0;
			for (int i = M+3; i < M+k+3; i++) knots(i) = 1;
			for(int i = 1;i < M-n+3; i++) {
				if(n > 1) {
					knots(n+i) = 0;
					for(int j = i-1; j < i+n-1; j++) knots(n+i) += bar_t(j);
					knots(n+i) /= n;
				} else knots(n+i) = bar_t(i);
			}
			return knots;

		} else {
			int n = k-1;
			Vector knots(M+1+k);
			for (int i = 0; i < n+1; i++) knots(i) = 0;
			for (int i = M+1; i < M+k+1; i++) knots(i) = 1;
			for(int i = 1;i < M-n+1; i++) {
				if(n > 1) {
					knots(n+i) = 0;
					for(int j = i; j < i+n; j++) knots(n+i) += bar_t(j);
					knots(n+i) /= n;
				} else knots(n+i) = bar_t(i);
			}
			return knots;
		}
	}

	namespace yml {

		void BSparambsp(const int &ndeg, const int &npoles, const int &npts, const Vector &bar_t, Vector &knots) {
			assert(npts == bar_t.size());

			knots.resize(npoles+ndeg+1);

			for(int i = 0; i < ndeg+1; i++) knots(i) = 0; //bar_t(0);

			Scalar d = (Scalar)npts / (npoles - ndeg);

			for (int i = 1; i < npoles-ndeg; i++) {
				int j = (int) (i * d);
				Scalar alpha = i * d - j;
				knots(ndeg+i) = bar_t(j-1) + alpha * (bar_t(j) - bar_t(j-1));
			}

			for (int i = npoles; i < npoles+ndeg+1; i++) knots(i) = 1.0; //bar_t(npts-1); 
		}
	}

	InterpolationBsplines::InterpolationBsplines() {}

	InterpolationBsplines::InterpolationBsplines(const PointSet &q, const int &k, const bool &natural_end) {
		int M = q.rows() - 1;

		Vector bar_t = centripetal_bar_t(q);
		adapt_bar_t(bar_t);
		if(verbose>=1) std::cerr << "bar_t = \n" << bar_t << std::endl;

		knots_ = average_knots(M, k, bar_t, natural_end);
		if(verbose>=1) std::cerr << "knots_ = \n" << knots_ << std::endl;

		if (natural_end) {

			Matrix equationL(M+3, M+3);
			for (int j = 0; j < M+3; j++) equationL(0, j) = d2B(j, k, bar_t(0));
			for (int j = 0; j < M+3; j++) equationL(1, j) = d2B(j, k, bar_t(M));
			for (int i = 2; i < M+3; i++) {
				for (int j = 0; j < M+3; j++) {
					equationL(i, j) = B(j, k, bar_t(i-2));
					//if(verbose>=2) std::cerr << "B(" << j << "," << k << "," << bar_t(i-2) << ") = " << equationL(i, j) << std::endl;  
				}
			}
			if(verbose>=1) std::cerr << "equationL = \n" << equationL << std::endl;

			Eigen::HouseholderQR<Matrix> qr(equationL);

			PointSet new_q(M+3, q.cols());
			new_q.topRows(2).setZero();
			new_q.block(2, 0, M+1, q.cols()) = q;

			p_ = qr.solve(new_q);

			if(verbose>=1) {
				std::cerr << "p_ = \n" << p_ << std::endl;
				std::cerr << "verify = \n" << equationL * p_ - new_q << std::endl; 
			}

		} else {

			Matrix equationL(M+1, M+1);
			for (int i = 0; i < M+1; i++) {
				for (int j = 0; j < M+1; j++) {
					equationL(i, j) = B(j, k, bar_t(i));
					//if(verbose>=2) std::cerr << "B(" << j << "," << k << "," << bar_t(i) << ") = " << equationL(i, j) << std::endl;  
				}
			}
			if(verbose>=1) std::cerr << "equationL = \n" << equationL << std::endl;

			Eigen::HouseholderQR<Matrix> qr(equationL);
			p_ = qr.solve(q);

			if(verbose>=1) std::cerr << "p_ = \n" << p_ << std::endl;
		}

		k_ = k;
	}

	Point InterpolationBsplines::eval(const Scalar &t) const {
		int M = p_.rows() - 1;

		Point rt = Point::Zero(p_.cols());
		for (int i = 0; i < M+1; i++) {
			Point temp = B(i, k_, t) * p_.row(i);
			rt += temp;
		}

		return rt;
	}

	ApproximationBsplines::ApproximationBsplines() {}

	ApproximationBsplines::ApproximationBsplines(const PointSet &q, const int &k, const int &npoles, const Scalar &lambda) {
		
		Vector bar_t = centripetal_bar_t(q);
		if(verbose>=1) std::cerr << "bar_t = \n" << bar_t << std::endl;

		int M = npoles - 2;

		yml::BSparambsp(k-1, M+2, bar_t.size(), bar_t, knots_);

		if(verbose>=1) std::cerr << "knots_ = \n" << knots_ << std::endl;
		
		int L = q.rows() - 2;

		Matrix N(L, M);
		for (int i = 0; i < L; i++) {
			for (int j = 0; j < M; j++) N(i, j) = B(j+1, k, bar_t(i+1));
		}
		if(verbose>=1) std::cerr << "N = \n" << N << std::endl;

		Matrix Q(L, q.cols());
		for (int i = 0; i < L; i++) Q.row(i) = q.row(i+1) - B(0, k, bar_t(i+1)) * q.row(0) - B(M+1, k, bar_t(i+1)) * q.row(L+1);
		if(verbose>=1) std::cerr << "Q = \n" << Q << std::endl;

		Matrix R = N.transpose() * Q;
		if(verbose>=1) std::cerr << "R = \n" << R << std::endl;

		p_.resize(M+2, q.cols());
		p_.row(0) = q.row(0);
		p_.row(M+1) = q.row(L+1);

		if (k > 2 && lambda > 0.0) {
			
			Matrix E(M, M);
			for (int i = 0; i < M; i++) {
				for (int j = i; j < M; j++) {
					yml::BScomptEij(k-1, M+2, knots_, 2, i+1, j+1, E(i,j));
					E(j, i) = E(i, j);
				}
			}
			if(verbose>=1) std::cerr << "E = \n" << E << std::endl;

			Matrix tilde_R(M, q.cols());
			tilde_R = R;
			for (int i = 0; i < k-1; i++) {
				Scalar tmp;
				yml::BScomptEij(k-1, M+2, knots_, 2, i+1, 0, tmp);
				tilde_R.row(i) = R.row(i) - lambda * p_.row(0) * tmp;
			}
			for (int i = M-k+1; i < M; i++) {
				Scalar tmp;
				yml::BScomptEij(k-1, M+2, knots_, 2, i+1, M+1, tmp);
				tilde_R.row(i) = R.row(i) - lambda * p_.row(M+1) * tmp; 
			}

			Matrix equationL(M, M);
			equationL = N.transpose() * N + lambda * E;

			Eigen::HouseholderQR<Matrix> qr(equationL);
			p_.middleRows(1, M) = qr.solve(tilde_R);
			if(verbose>=1) std::cerr << "p_ = \n" << p_ << std::endl;

		} else {

			Matrix equationL(M, M);
			equationL = N.transpose() * N;

			Eigen::HouseholderQR<Matrix> qr(equationL);
			p_.middleRows(1, M) = qr.solve(R);
			if(verbose>=1) std::cerr << "p_ = \n" << p_ << std::endl;
		}

		k_ = k;
	}

	Point ApproximationBsplines::eval(const Scalar &t) const {
		int M = p_.rows() - 1;

		Point rt = Point::Zero(p_.cols());
		for (int i = 0; i < M+1; i++) {
			Point temp = B(i, k_, t) * p_.row(i);
			rt += temp;
		}

		return rt;
	}

	namespace yml {
		void BScomptEij(const int &ndeg, const int &npoles, const Vector &knots, const int &nder, const int &ii, const int &jj, Scalar &dValue) {
			assert(ndeg <= 7);

			BSplines bsplines(knots);
			BSplinesFunction dxB_i_k(bsplines, ii, ndeg+1, nder);
			BSplinesFunction dxB_j_k(bsplines, jj, ndeg+1, nder);
			TensorBSplinesFunction tensorBB(dxB_i_k, dxB_j_k);
			dValue = 0.0;
 			for(int l = std::max(ndeg+1, std::max(ii+1, jj+1)); l <= (ii+ndeg+1) && l <= (jj+ndeg+1); l++) dValue += gaussion_quadrature( tensorBB, ndeg - nder + 1, knots(l-1), knots(l));
		}
	}

	BSplinesFunction::BSplinesFunction(const BSplines &bsplines, const int &i, const int &k, const int &nder) : bsplines_(bsplines), i_(i), k_(k), nder_(nder) {}

	Scalar BSplinesFunction::operator()(const Scalar &t) const {
		assert(nder_ <= 2);
		switch(nder_) {
		case 0: return bsplines_.B(i_, k_, t);
		case 1: return bsplines_.dB(i_, k_, t);
		case 2: return bsplines_.d2B(i_, k_, t);
		}
		return 0.0;
	}

	TensorBSplinesFunction::TensorBSplinesFunction(const BSplinesFunction &bspfunc1, const BSplinesFunction &bspfunc2) : bspfunc1_(bspfunc1), bspfunc2_(bspfunc2) {}

	Scalar TensorBSplinesFunction::operator()(const Scalar &t) const {
		return bspfunc1_(t) * bspfunc2_(t);
	}

	GFunction::GFunction(const Function &f, const Scalar &a, const Scalar &b) : f_(f), a_(a), b_(b) {}

	Scalar GFunction::operator()(const Scalar &t) const {
		return f_( ( (a_+b_) + (b_-a_) * t) * 0.5 );
	}

	Scalar gaussion_quadrature(const Function &f, const int &n, const Scalar &a, const Scalar &b) {

		static Scalar x1[1] = {0.0}, x2[2] = {-sqrt(3.0)/3, sqrt(3.0)/3}, x3[3] = {-sqrt(3.0/5), 0, sqrt(3.0/5)},
			x4[4] = {-sqrt( (15 + 2 * sqrt(30.0))/35 ), -sqrt( (15 - 2 * sqrt(30.0))/35 ), sqrt( (15 - 2 * sqrt(30.0))/35 ), sqrt( (15 + 2 * sqrt(30.0))/35 ) };
		static Scalar a1[1] = {2.0}, a2[2] = {1.0, 1.0}, a3[3] = {5.0/9, 8.0/9, 5.0/9},
			a4[4] = { (18 - sqrt(30.0))/36, (18 + sqrt(30.0))/36, (18 + sqrt(30.0))/36, (18 - sqrt(30.0))/36};

		assert(n <= 4);

		GFunction g(f, a, b);
		Scalar rt;
		switch(n) {
		case 1: {
				rt = g(x1[0]) * a1[0];
				break;
			}
		case 2: {
				rt = g(x2[0]) * a2[0] + g(x2[1]) * a2[1];
			}
		case 3: {
				rt = g(x3[0]) * a3[0] + g(x3[1]) * a3[1] + g(x3[2]) * a3[2];
				break;
			}
		case 4: {
				rt = g(x4[0]) * a4[0] + g(x4[1]) * a4[1] + g(x4[2]) * a4[2] + g(x4[3]) * a4[3];
				break;
			}
		}
		return rt * (b - a) * 0.5;
	}


	namespace yml {

		void BScompnmat(const int &npts, const Vector &bar_t, const int &ndeg, const int &npoles, const Vector &knots, Matrix &Nmat) {

			BSplines bsplines(knots);
			
			Nmat.resize(npts-2, npoles-2);
			for (int i = 0; i < npts-2; i++) {
				for (int j = 0; j < npoles-2; j++) {
					Nmat(i, j) = bsplines.B(j+1, ndeg+1, bar_t(i+1));
				}
			}
		}

		void BSprodnmat(const int &npts, const int &npoles, const Matrix &Nmat, Matrix &Amat) {

			Amat = Nmat.transpose() * Nmat;
		}

		void BScomprvec(const int &npts, const int &ndim, const PointSet &Qpts, const Vector &bar_t, const Matrix &Nmat, const int &ndeg, const int &npoles, const Vector &knots, PointSet &Rvec) {
			
			BSplines bsplines(knots);

			PointSet Q(npts-2, ndim);
			for (int i = 0; i < npts-2; i++) {
				Q.row(i) = Qpts.row(i+1) - bsplines.B(0, ndeg+1, bar_t(i+1)) * Qpts.row(0) - bsplines.B(npoles-1, ndeg+1, bar_t(i+1)) * Qpts.row(npts-1);
			}

			Rvec = Nmat.transpose() * Q;
		}

		void BSbndcholmt(const int &npoleu, const int &dim, const Matrix &Amat, const PointSet &Rvec, PointSet &sol) {
			
			Eigen::HouseholderQR<Matrix> qr(Amat);
			sol = qr.solve(Rvec);
		}

		void BSlstfitptsf(const int &nptu, const int &nptv, const PointSet &Qpts, const Vector &bar_u, const Vector &bar_v, const int &ndegu, const int &ndegv, const int &npoleu, const int &npolev,
			PointSet &poles, Vector &uknots, Vector &vknots) {
				
				int dim = Qpts.cols();

				int npts_max = std::max(nptu, nptv);
				int npole_max = std::max(npoleu, npolev);
				int ndeg_max = std::max(ndegu, ndegv);
				Matrix Nmat(npts_max, ndeg_max);
				Matrix Amat(npole_max, npole_max);
				PointSet Rvec(npole_max, dim);
				PointSet poles_temp(nptv * npoleu, dim);

				BSparambsp(ndegu, npoleu, nptu, bar_u, uknots);
				BSparambsp(ndegv, npolev, nptv, bar_v, vknots);

				BScompnmat(nptu, bar_u, ndegu, npoleu, uknots, Nmat);
				BSprodnmat(nptu, npoleu, Nmat, Amat);
				for (int j = 0; j < nptv; j++) {
					BScomprvec(nptu, dim, Qpts.middleRows(j * nptu, nptu), bar_u, Nmat, ndegu, npoleu, uknots, Rvec);
					PointSet sol(npoleu-2, dim);
					BSbndcholmt(npoleu-2, dim, Amat, Rvec, sol);

					poles_temp.row(j) = Qpts.middleRows(j * nptu, nptu).row(0);
					for(int i = 0; i < npoleu-2; i++) poles_temp.row( nptv * (i+1) + j) = sol.row(i);
					poles_temp.row(nptv * (npoleu - 1) + j) = Qpts.middleRows(j * nptu, nptu).row(nptu-1);
				}

				poles.resize(npoleu * npolev, dim);

				BScompnmat(nptv, bar_v, ndegv, npolev, vknots, Nmat);
				BSprodnmat(nptv, npolev, Nmat, Amat);
				for (int i = 0; i < npoleu; i++) {
					BScomprvec(nptv, dim, poles_temp.middleRows(i * nptv, nptv), bar_v, Nmat, ndegv, npolev, vknots, Rvec);
					PointSet sol(npolev-2, dim);
					BSbndcholmt(npolev-2, dim, Amat, Rvec, sol);

					poles.row(i) = poles_temp.middleRows(i * nptv, nptv).row(0);
					for(int j = 0; j < npolev-2; j++) poles.row(npoleu * (j+1) + i) = sol.row(j);
					poles.row( npoleu * (npolev-1) + i ) = poles_temp.middleRows(i * nptv, nptv).row(nptv-1);
				}
		}
	}

	TensorBSplines::TensorBSplines() {}

	TensorBSplines::TensorBSplines(const int &nptu, const int &nptv, const PointSet &Qpts, const Vector &bar_u, const Vector &bar_v, const int &ndegu, const int &ndegv, const int &npoleu, const int &npolev) {
		yml::BSlstfitptsf(nptu, nptv, Qpts, bar_u, bar_v, ndegu, ndegv, npoleu, npolev, poles_, uknots_, vknots_);
		npoleu_ = npoleu;
		npolev_ = npolev;
		ku_ = ndegu+1;
		kv_ = ndegv+1;
	}

	Point TensorBSplines::eval(const Scalar &u, const Scalar &v) const {
		
		int dim = poles_.cols();

		BSplines bsplines_u(uknots_);
		BSplines bsplines_v(vknots_);
		
		Point rt = Point::Zero(dim);
		for (int j = 0; j < npolev_; j++) {
			for (int i = 0; i < npoleu_; i++) {
				rt += poles_.row(j * npoleu_ + i) * bsplines_u.B(i, ku_, u) * bsplines_v.B(j, kv_, v);
			}
		}

		return rt;
	}

	TensorBSplines::TensorBSplines(const int &npts, const PointSet &Qpts, const Vector &bar_u, const Vector &bar_v, const int &ndegu, const int &ndegv, 
		const int &npoleu, const int &npolev, const Scalar &lambda) {

		Vector temp_bar_u(bar_u);
		std::sort(temp_bar_u.data(), temp_bar_u.data() + temp_bar_u.size());

		Vector temp_bar_v(bar_v);
		std::sort(temp_bar_v.data(), temp_bar_v.data() + temp_bar_v.size());

		yml::BSparambsp(ndegu, npoleu, temp_bar_u.size(), temp_bar_u, uknots_);
		yml::BSparambsp(ndegv, npolev, temp_bar_v.size(), temp_bar_v, vknots_);

		if(verbose>=1) {
			std::cerr << "uknots_ = \n" << uknots_ << std::endl;
			std::cerr << "vknots_ = \n" << vknots_ << std::endl;
		}

		BSplines bsplines_u(uknots_);
		BSplines bsplines_v(vknots_);

		Matrix Bmat(npts, npoleu * npolev);

		for (int l = 0; l < npts; l++) {
			for (int k = 0; k < npoleu * npolev; k++) {
				int i = k % npoleu, j = k / npoleu;
				Bmat(l, k) = bsplines_u.B(i, ndegu, bar_u(l)) * bsplines_v.B(j, ndegv, bar_v(l));
			}
		}

		if(verbose>=1) std::cerr << "Bmat = \n" << Bmat << std::endl;

		if ( (ndegu+ndegv > 1) && (lambda > 0.0) ) {

			Matrix Emat(npoleu * npolev, npoleu * npolev), Dmat(npoleu * npolev, npoleu * npolev), Gmat(npoleu * npolev, npoleu * npolev), Fmat(npoleu * npolev, npoleu * npolev);

			Matrix NNmat_u_0, NNmat_u_1, NNmat_u_2;
			yml::BScomptnnmat(ndegu, npoleu, uknots_, 0, NNmat_u_0);
			yml::BScomptnnmat(ndegu, npoleu, uknots_, 1, NNmat_u_1);
			yml::BScomptnnmat(ndegu, npoleu, uknots_, 2, NNmat_u_2);

			Matrix NNmat_v_0, NNmat_v_1, NNmat_v_2;
			yml::BScomptnnmat(ndegv, npolev, vknots_, 0, NNmat_v_0);
			yml::BScomptnnmat(ndegv, npolev, vknots_, 1, NNmat_v_1);
			yml::BScomptnnmat(ndegv, npolev, vknots_, 2, NNmat_v_2);

			for (int kj = 0; kj < npoleu * npolev; kj++) {
				for (int ki = 0; ki < npoleu * npolev; ki++) {
					int i = kj % npoleu, j = kj / npoleu;
					int r = ki % npoleu, s = ki / npoleu;

					Dmat(kj, ki) = NNmat_u_2(i, r) * NNmat_v_0(j, s);
					Gmat(kj, ki) = NNmat_u_1(i, r) * NNmat_v_1(j, s);
					Fmat(kj, ki) = NNmat_u_0(i, r) * NNmat_v_2(j, s);
				}
			}

			Emat = Dmat + 2 * Gmat + Fmat;

			Eigen::HouseholderQR<Matrix> qr(Bmat.transpose() * Bmat + lambda * Emat);
			poles_ = qr.solve(Bmat.transpose() * Qpts);
		
		} else {
			Eigen::HouseholderQR<Matrix> qr(Bmat.transpose() * Bmat);
			poles_ = qr.solve(Bmat.transpose() * Qpts);
		}

		if(verbose>=1) std::cerr << "poles_ = \n" << poles_ << std::endl;

		npoleu_ = npoleu;
		npolev_ = npolev;

		ku_ = ndegu+1;
		kv_ = ndegv+1;
	}

	namespace yml {

		void BScomptnnmat(const int &ndeg, const int &npoles, const Vector &knots, const int &nder, Matrix &NNmat) {
			
			NNmat = Matrix::Zero(npoles, npoles);

			for (int j = 0; j < npoles; j++) {
				int n = j * npoles;
				int end = std::min(j+ndeg+1, npoles);
				for (int s = std::max(0, j-ndeg); s < end; s++) {
					int m = s * npoles;
					if (s < j) {
						NNmat(s + n) = NNmat(m + j);
					}
					else {
						BScomptEij(ndeg, npoles, knots, nder, j, s, NNmat(s+n));
					}
				}
			}
		}
	}
	


}