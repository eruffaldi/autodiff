/**
 * Design Choice 1: Eigen/Dense containing sym objects 
 * Design Choice 2: interface similar to Eigen with special markers
 *
 * Can we support both?
 */
#pragma once
#include "sympp.hpp"
#include <Eigen/Dense>

#if 0
namespace Eigen
{
template<> struct NumTraits<sym> : GenericNumTraits<sym>
{
    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 10,
      AddCost = 10,
      MulCost = 40
    };

    typedef sym Real;
    typedef sym NonInteger;
    typedef sym Nested;
    
    inline static Real Pi       (long Precision = 0)     {    return pi();        }
    inline static Real Euler    (long Precision = 0)     {    return e();    }
    inline static Real Log2     (long Precision = 0)     {    return log2();     }
};
#endif

enum class MatrixContent {
	Zero,
	Identity,
	Generic 
};

enum class MatrixStorage { 
	Scalar, // diagonal with all same elements
	Uniform, // all same values (equivalent to Scalar if scalar size)
	Diagonal, // diagonal
	UpperTriangular, 
	LowerTriangular, 
	Symmetric, 
	Generic 
};

/// default is scalar
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
	virtual ~imsym() {}
	//virtual std::shared_ptr<isym> diff(int index) const = 0;
	virtual std::string sig() const = 0;
	virtual void print(std::ostream & os) const = 0;
	virtual Eigen::MatrixXd operator()() const = 0;
	virtual std::shared_ptr<imsym> diff(int p) const = 0;
	virtual int nparents() const = 0;
	virtual std::shared_ptr<imsym> parent(int p) const = 0;
	virtual bool isconst() const { return false; }
	virtual bool set(Eigen::MatrixXd) { return false; }
	virtual msym parentadjointS(int i) const = 0;
	virtual Eigen::MatrixXd  parentadjointN(int i) const = 0;

	imsym tosym() const;
};


/// valued type
struct msym
{
	/// variable
	explicit msym(std::string name, int rows, int cols, MatrixStorage s = MatrixStorage::Scalar,MatrixContent cc = MatrixContent::Generic );

	/// variable
	explicit msym(std::string name, matrixspec spec);

	/// valued matrix
	explicit msym(const Eigen::MatrixXd & m);

	void print(std::ostream & os) const { p_->print(os); }

#if 0
	msym &operator += (const msym & a) ;
	msym &operator *= (const msym & a) ;
	msym &operator -= (const msym & a) ;
	msym operator - () const;
#endif

	/// evaluate
	Eigen::MatrixXd val() const { return (*p_)(); }

	int nparents() const { return p_->nparents(); }
	msym parent(int i) const { return (p_->parent(i))->tosym(); }
	msym diff(int iparent) const { return (p_->diff(iparent))->tosym(); }

	//sym solveconst();

	bool operator== (const msym &x) const { return x.p_.get() == p_.get(); }

private:
	std::shared_ptr<imsym> p_;
};


inline std::ostream & operator << (std::ostream & os, const msym & x) 
{
	x.print(os);
	return os;
}

#if 0
msym operator + (const msym & a, const msym & b);
msym operator - (const msym & a, const msym & b);
msym operator * (const msym & a, const msym & b);

msym operator + (double a, const msym & b);
msym operator - (double a, const msym & b);
msym operator / (double a, const msym & b);
msym operator * (double a, const msym & b);
msym operator + (const msym & a, double b);
msym operator - (const msym & a, double b);
msym operator / (const msym & a, double b);
msym operator * (const msym & a, double b);

msym operator + (sym a, const msym & b);
msym operator - (sym a, const msym & b);
msym operator / (sym a, const msym & b);
msym operator * (sym a, const msym & b);
msym operator + (const msym & a, sym b);
msym operator - (const msym & a, sym b);
msym operator / (const msym & a, sym b);
msym operator * (const msym & a, sym b);

msym transpose(msym x);
//msym trace(msym x);
msym frobenius(msym x);
//msym operator / (const msym & a, const msym & b);
msym det(msym x);
msym invert(msym x);
#endif

#if 0
	/// valued variable
	sym(std::string name, const Eigen::MatrixXd & value);
/// constant
sym::sym(const Eigen::MatrixXd & value): p_(std::make_shared<vconstmat>(value))
{

}

#endif