/**
 * C++ Automatic Diff exercise: reverse mode, Matrix, dynamic
 * Emanuele Ruffaldi
 *
 * Fundamentals from stevenrennie+Olen IBM 
 * - Slides: http://www.stevenrennie.com/papers/steve_rennie_matdiff_7_19_2012.pdf
 * - Paper: http://www.stevenrennie.com/papers/jfd.pdf
 */
#pragma once
#include <memory>
#include <string>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

struct msym;

/// default is scalar (NOT USED SO-FAR)
struct matrixspec
{
	matrixspec() {}

	matrixspec(int r,int c,MatrixStorage s = MatrixStorage::Scalar,MatrixContent cc = MatrixContent::Generic): rows(r),cols(c),storage(s), content(cc) {}
	MatrixStorage storage = MatrixStorage::Scalar;
	MatrixContent content = MatrixContent::Generic;
	int rows = 1;
	int cols = 1;

	bool isscalar() const { return rows == 1 && cols == 1; }

	bool isvector() const { return rows == 1 || cols == 1; }

	bool issquare() const { return rows == cols; }

	bool isidentity() const { return issquare() && content == MatrixContent::Identity; }

	matrixspec transpose() const { return matrixspec(cols,rows,storage,content); }
};

/// base class
struct imsym : public std::enable_shared_from_this<imsym>
{
	using mat_t = Eigen::MatrixXd;

	virtual ~imsym() {}
	//virtual std::shared_ptr<imsym> diff(int index) const = 0;
	virtual std::string sig() const = 0;
	virtual void print(std::ostream & os) const = 0;
	virtual mat_t operator()() const = 0;
	virtual bool set(mat_t v)  = 0;
	//virtual std::shared_ptr<imsym> diff(int p) const = 0;
	virtual int nparents() const = 0;
	virtual std::shared_ptr<imsym> parent(int p) const = 0;
	virtual bool isconst() const { return false; }

	/// returns the symbolic adjoint contribution to the i-th parent. 
	/// NOTE: in the case of matrix expression graph the adjoint is not simply: diff * adjoint but a more general expression involving the adjoint	
//	virtual msym parentadjointS(int i, msym adjoint) const = 0;

	/// returns the numeric adjoint
	virtual mat_t  parentadjointN(int i, const mat_m & adjoint) const = 0;


	msym tomsym() const;
};

/// valued type
struct msym
{
	using mat_t = typename imsym::mat_t;
	/// wrap
	msym(const std::shared_ptr<imsym> & pp): p_(pp) {}

	/// variable
	explicit msym(std::string name);

	/// valued variable
	msym(std::string name, mat_t d);

	/// constant
	explicit msym(double d);

	/// constant
	explicit msym(mat_t d);

	/// constant
	explicit msym(int d);

	virtual bool set(mat_t v) {return p_->set(v); }

	void print(std::ostream & os) const { p_->print(os); }

	msym &operator += (const msym & a) ;
	msym &operator *= (const msym & a) ;
	msym &operator /= (const msym & a) ;
	msym &operator -= (const msym & a) ;
	msym operator - () const;

	/// evaluate
	mat_t val() const { return (*p_)(); }

	/// related terms
	int nparents() const { return p_->nparents(); }
	msym parent(int i) const { return (p_->parent(i))->tomsym(); }
//	msym diff(int iparent) const { return (p_->diff(iparent))->tomsym(); }

	/// returns the numeric adjoint
	virtual mat_t  parentadjointN(int i, const mat_m & adjoint) const { return p_->parentadjointN(i,adjoint); }

	//msym solveconst();

	std::shared_ptr<imsym> p_;	

	bool operator== (const msym &x) const { return x.p_.get() == p_.get(); }
};

inline std::ostream & operator << (std::ostream & os, const msym & x) 
{
	x.print(os);
	return os;
}

msym operator + (const msym & a, const msym & b);
msym operator - (const msym & a, const msym & b);
msym operator / (const msym & a, const msym & b);
msym operator * (const msym & a, const msym & b);

msym operator + (double a, const msym & b);
msym operator - (double a, const msym & b);
msym operator / (double a, const msym & b);
msym operator * (double a, const msym & b);

msym operator + (const msym & a, double b);
msym operator - (const msym & a, double b);
msym operator / (const msym & a, double b);
msym operator * (const msym & a, double b);

msym operator + (msym::mat_t a, const msym & b);
msym operator - (msym::mat_t a, const msym & b);
msym operator * (msym::mat_t a, const msym & b);

msym operator + (const msym & a, msym::mat_t b);
msym operator - (const msym & a, msym::mat_t b);
msym operator * (const msym & a, msym::mat_t b);

// external unique operations
msym one();
msym negone();
msym zero();
msym two();
msym half();
msym msym_e();
msym msym_pi();

template <class T>
class pi
{
};

template <>
class pi<msym>
{
public:
	operator msym() { return msym_pi(); }
};
template <>
class pi<double>
{
public:
	operator double() { return M_PI; }
};

msym log2();
msym sqrt(msym);
msym log(msym);
msym cos(msym);
msym exp(msym);
msym sin(msym);
msym pow(msym x, double q) ;
msym pow(msym x1, msym x2);
msym invert(msym x);

#if 0
/// msymbolic jacobian
struct jacobmsym
{

	jacobmsym(msym root,std::initializer_list<msym> vars);

	/// number of gradients
	int gradients() const { return targets_.size(); } 

	/// i-th target
	msym gradientVar(int i) const { return targets_[i]; }

	/// i-th
	msym gradient(int i) const { return gradient(targets_[i]); }

	/// by msymbol
	msym gradient(msym s) const;


private:
	std::vector<msym> targets_;
	std::unordered_map<std::shared_ptr<imsym>, msym > adjoint_;
	void descend(msym v,bool isroot);
};
#endif

/// numeric jacobian
struct jacobnum
{

	jacobnum(msym root,std::initializer_list<msym> vars);

	void update(); // when new values are set to the depending vars
	int gradients() const { return targets_.size(); } 
	mat_t gradient(int i) const { return gradient(targets_[i]); }
	msym gradientVar(int i) const { return targets_[i]; }
	mat_t gradient(msym s) const ;
private:
	void updateTargetSize(); 
	Eigen::Vector2i inputsize;
	
	msym root_;
	std::vector<msym> targets_;
	std::unordered_map<std::shared_ptr<imsym>, mat_t > adjoint_;
	std::unordered_map<std::shared_ptr<imsym>, std::vector<msym> > diff_; 
	using pti = typename std::unordered_map<std::shared_ptr<imsym>, mat_t >::iterator;

	/// recursive descent of the expression tree down to leaves for computing the derivative expression wrt each term
	/// the derivative expression is stored in diff
	/// In dx/dy = dx/dw dw/dy this builds dx/dw
	//void descendG(msym v);

	// NOTE: after FIRST run NO more allocations needed
	//void descendR(msym v,pti it);
};





//msym tan(msym);
//msym atan(msym);
//msym atan2(msym,msym); // TODO: mixed variant