#pragma once
#include <memory>
#include <string>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

struct sym;

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
struct isym : public std::enable_shared_from_this<isym>
{
	virtual ~isym() {}
	//virtual std::shared_ptr<isym> diff(int index) const = 0;
	virtual std::string sig() const = 0;
	virtual void print(std::ostream & os) const = 0;
	virtual double operator()() const = 0;
	virtual std::shared_ptr<isym> diff(int p) const = 0;
	virtual int nparents() const = 0;
	virtual std::shared_ptr<isym> parent(int p) const = 0;
	virtual bool isconst() const { return false; }

	sym tosym() const;

	matrixspec spec; /// default is scalar
};

/// valued type
struct sym
{
	/// wrap
	sym(const std::shared_ptr<isym> & pp): p_(pp) {}

	/// variable
	explicit sym(std::string name);

	/// valued variable
	sym(std::string name, const double d);

	/// valued variable
	sym(std::string name, const Eigen::MatrixXd & value);

	/// constant
	explicit sym(double d);

	/// constant
	explicit sym(int d);

	/// constant
	explicit sym(const Eigen::MatrixXd & value);

	void print(std::ostream & os) const { p_->print(os); }

	sym &operator += (const sym & a) ;
	sym &operator *= (const sym & a) ;
	sym &operator /= (const sym & a) ;
	sym &operator -= (const sym & a) ;
	sym operator - () const;

	/// evaluate
	double val() const { return (*p_)(); }

	int nparents() const { return p_->nparents(); }
	sym parent(int i) const { return (p_->parent(i))->tosym(); }
	sym diff(int iparent) const { return (p_->diff(iparent))->tosym(); }

	matrixspec spec() const { return p_->spec; }
	sym solveconst();

	std::shared_ptr<isym> p_;	

	bool operator== (const sym &x) const { return x.p_.get() == p_.get(); }
};

inline std::ostream & operator << (std::ostream & os, const sym & x) 
{
	x.print(os);
	return os;
}


sym operator + (const sym & a, const sym & b);
sym operator - (const sym & a, const sym & b);
sym operator / (const sym & a, const sym & b);
sym operator * (const sym & a, const sym & b);

sym operator + (double a, const sym & b);
sym operator - (double a, const sym & b);
sym operator / (double a, const sym & b);
sym operator * (double a, const sym & b);

sym operator + (const sym & a, double b);
sym operator - (const sym & a, double b);
sym operator / (const sym & a, double b);
sym operator * (const sym & a, double b);

// external unique operations
sym one();
sym negone();
sym zero();
sym two();
sym e();
sym pi();
sym sqrt(sym);
sym log(sym);
sym cos(sym);
sym exp(sym);
sym sin(sym);
sym pow(sym x, double q) ;
sym pow(sym x1, sym x2);
sym invert(sym x);

/// matrix specific
sym transpose(sym x);
sym det(sym x);
sym trace(sym x);
sym frobenius(sym x);

/// symbolic jacobian
struct jacob
{

	jacob(sym root,std::initializer_list<sym> vars);

	int gradients() const { return targets_.size(); } 
	sym gradient(int i) const { return gradient(targets_[i]); }
	sym gradient(sym s) const;
	sym gradientVar(int i) const { return targets_[i]; }

private:
	std::vector<sym> targets_;
	std::unordered_map<std::shared_ptr<isym>, sym > adjoint_;
	void descend(sym v,bool isroot);
};

/// symbolic jacobian
struct jacobnum
{

	jacobnum(sym root,std::initializer_list<sym> vars);

	void update(); // when new values are set to the depending vars
	int gradients() const { return targets_.size(); } 
	double gradient(int i) const { return gradient(targets_[i]); }
	sym gradientVar(int i) const { return targets_[i]; }
	double gradient(sym s) const ; 

private:
	sym root_;
	std::vector<sym> targets_;
	std::unordered_map<std::shared_ptr<isym>, double > adjoint_;
	std::unordered_map<std::shared_ptr<isym>, std::vector<sym> > diff_;
	void descendG(sym v);
	void descendR(sym v);
};


//sym tan(sym);
//sym atan(sym);
//sym atan2(sym,sym); // TODO: mixed variant