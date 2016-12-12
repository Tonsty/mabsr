#ifndef MABSR_BSPLINES_H
#define MABSR_BSPLINES_H

#ifndef LIBBSPLINES_API
#ifdef _WIN32
#ifdef LIBBSPLINES_DYNAMIC
#if LIBBSPLINES_BUILD
#define LIBBSPLINES_API __declspec(dllexport)
#else 
#define LIBBSPLINES_API __declspec(dllimport)
#endif
#else
#define LIBBSPLINES_API
#endif
#else
#define LIBBSPLINES_API
#endif
#endif

#include <Eigen/Dense>

namespace mabsr {

	typedef float Scalar;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
	typedef Vector Point;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> PointSet;
	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

	class BSplines {
	public:
		LIBBSPLINES_API BSplines();

		LIBBSPLINES_API BSplines(const Vector &_knots);

		LIBBSPLINES_API Scalar B(const int &i, const int &k, const Scalar &t) const;

		LIBBSPLINES_API Scalar dB(const int &i, const int &k, const Scalar &t) const;

		LIBBSPLINES_API Scalar d2B(const int &i, const int &k, const Scalar &t) const;

		LIBBSPLINES_API Scalar w(const int &i, const int &k, const Scalar &t) const;

		LIBBSPLINES_API Scalar X(const int &i, const Scalar &t) const;

	protected:

		Vector knots_;
	};

	Vector equally_spaced_bar_t(const PointSet &q);

	Vector chordal_bar_t(const PointSet &q);

	Vector centripetal_bar_t(const PointSet &q);

	LIBBSPLINES_API void adapt_bar_t(Vector &bar_t, const Scalar &scale = 0.999999, const Scalar &shift = 0.0000005);

	Vector equally_spaced_knots(const int &M, const int &k, const bool &natural_end = false);

	Vector average_knots(const int &M, const int &k, const Vector &bar_t, const bool &natural_end = false);

	namespace yml {
		LIBBSPLINES_API void BSparambsp(const int &ndeg, const int &npoles, const int &npts, const Vector &bar_t, Vector &knots);
	}


	class InterpolationBsplines: public BSplines {
	public:
		LIBBSPLINES_API InterpolationBsplines();

		LIBBSPLINES_API InterpolationBsplines(const PointSet &q, const int &k, const bool &_natural_end = false);

		LIBBSPLINES_API Point eval(const Scalar &t) const;

	protected:

		PointSet p_;
		int k_;
	};

	class ApproximationBsplines: public BSplines {
	public:
		LIBBSPLINES_API ApproximationBsplines();

		LIBBSPLINES_API ApproximationBsplines(const PointSet &q, const int &k, const int &M, const Scalar &lambda = 0.0);

		LIBBSPLINES_API Point eval(const Scalar &t) const;

	protected:
		PointSet p_;
		int k_;
	};

	namespace yml {
		void BScomptEij(const int &ndeg, const int &npoles, const Vector &knots, const int &nder, const int &ii, const int &jj, Scalar &dValue);
	}

	class Function {
	public:
		virtual Scalar operator()(const Scalar &t) const = 0;
	};

	class BSplinesFunction : public Function{
	public:
		BSplinesFunction(const BSplines &bsplines, const int &i, const int &k, const int &nder = 0);
		virtual Scalar operator()(const Scalar &t) const;

	protected:
		const BSplines &bsplines_;
		const int i_;
		const int k_;
		const int nder_;
	};

	class TensorBSplinesFunction : public Function{
	public:
		TensorBSplinesFunction(const BSplinesFunction &bspfunc1, const BSplinesFunction &bspfunc2);
		virtual Scalar operator()(const Scalar &t) const;

	protected:
		const BSplinesFunction &bspfunc1_;
		const BSplinesFunction &bspfunc2_;
	};

	class GFunction : Function{
	public:
		GFunction(const Function &f, const Scalar &a, const Scalar &b);
		virtual Scalar operator()(const Scalar &t) const;

	protected:
		const Function &f_;
		const Scalar &a_;
		const Scalar &b_;
	};

	Scalar gaussion_quadrature(const Function &f, const int &n = 3, const Scalar &a = -1.0, const Scalar &b = 1.0);

	namespace yml {

		void BScompnmat(const int &npts, const Vector &bar_t, const int &ndeg, const int &npoles, const Vector &knots, Matrix &Nmat);

		void BSprodnmat(const int &npts, const int &npoles, const Matrix &Nmat, Matrix &Amat);

		void BScomprvec(const int &npts, const int &ndim, const PointSet &Qpts, const Vector &bar_t, const Matrix &Nmat, const int &ndeg, const int &npoles, const Vector &knots, PointSet &Rvec);

		void BSbndcholmt(const int &npoleu, const int &dim, const Matrix &Amat, const PointSet &Rvec, PointSet &sol);

		void BSlstfitptsf(const int &nptu, const int &nptv, const PointSet &Qpts, const Vector &bar_u, const Vector &bar_v, const int &ndegu, const int &ndegv, const int &npoleu, const int &npolev,
			PointSet &poles, Vector &uknots, Vector &vknots);
	};

	class TensorBSplines {
	public:
		LIBBSPLINES_API TensorBSplines();
		
		LIBBSPLINES_API TensorBSplines(const int &nptu, const int &nptv, const PointSet &Qpts, const Vector &bar_u, const Vector &bar_v, const int &ndegu, const int &ndegv, const int &npoleu, const int &npolev);

		LIBBSPLINES_API TensorBSplines(const int &npts, const PointSet &Qpts, const Vector &bar_u, const Vector &bar_v, const int &ndegu, const int &ndegv, const int &npoleu, const int &npolev, const Scalar &lambda = 0.0);

		LIBBSPLINES_API Point eval(const Scalar &u, const Scalar &v) const;

	protected:
		PointSet poles_;
		int npoleu_;
		int npolev_;
		int ku_;
		int kv_;
		Vector uknots_;
		Vector vknots_;
	};

	namespace yml {

		void BScomptnnmat(const int &ndeg, const int &npoles, const Vector &knots, const int &nder, Matrix &NNmat);
	}

};

#endif
