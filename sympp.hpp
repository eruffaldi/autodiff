/**
 * C++ Automatic Diff exercise: reverse mode, scalar only
 * Emanuele Ruffaldi
 */
#pragma once
#include <memory>
#include <string>
#include <unordered_map>
#include <iostream>
#include <vector>

struct sym;

/// base class
struct isym : public std::enable_shared_from_this<isym>
{
	virtual ~isym() {}
	//virtual std::shared_ptr<isym> diff(int index) const = 0;
	virtual std::string sig() const = 0;
	virtual void print(std::ostream & os) const = 0;
	virtual double operator()() const = 0;
	virtual bool set(double v)  = 0;
	virtual std::shared_ptr<isym> diff(int p) const = 0;
	virtual int nparents() const = 0;
	virtual std::shared_ptr<isym> parent(int p) const = 0;
	virtual bool isconst() const { return false; }

	sym tosym() const;
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

	/// constant
	explicit sym(double d);

	/// constant
	explicit sym(int d);

	virtual bool set(double v) {return p_->set(v); }

	void print(std::ostream & os) const { p_->print(os); }

	sym &operator += (const sym & a) ;
	sym &operator *= (const sym & a) ;
	sym &operator /= (const sym & a) ;
	sym &operator -= (const sym & a) ;
	sym operator - () const;

	/// evaluate
	double val() const { return (*p_)(); }

	/// related terms
	int nparents() const { return p_->nparents(); }
	sym parent(int i) const { return (p_->parent(i))->tosym(); }
	sym diff(int iparent) const { return (p_->diff(iparent))->tosym(); }

	//sym solveconst();

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
sym half();
sym sym_e();
sym sym_pi();

template <class T>
class pi
{
};

template <>
class pi<sym>
{
public:
	operator sym() { return sym_pi(); }
};
template <>
class pi<double>
{
public:
	operator double() { return M_PI; }
};

sym log2();
sym sqrt(sym);
sym log(sym);
sym cos(sym);
sym exp(sym);
sym sin(sym);
sym pow(sym x, double q) ;
sym pow(sym x1, sym x2);
sym invert(sym x);


/// symbolic jacobian
struct jacobsym
{

	jacobsym(sym root,std::initializer_list<sym> vars);

	/// number of gradients
	int gradients() const { return targets_.size(); } 

	/// i-th target
	sym gradientVar(int i) const { return targets_[i]; }

	/// i-th
	sym gradient(int i) const { return gradient(targets_[i]); }

	/// by symbol
	sym gradient(sym s) const;


private:
	std::vector<sym> targets_;
	std::unordered_map<std::shared_ptr<isym>, sym > adjoint_;
	void descend(sym v,bool isroot);
};

/// numeric jacobian
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
	using pti = typename std::unordered_map<std::shared_ptr<isym>, double >::iterator;

	/// recursive descent of the expression tree down to leaves for computing the derivative expression wrt each term
	/// the derivative expression is stored in diff
	/// In dx/dy = dx/dw dw/dy this builds dx/dw
	void descendG(sym v);

	// NOTE: after FIRST run NO more allocations needed
	void descendR(sym v,pti it);
};





//sym tan(sym);
//sym atan(sym);
//sym atan2(sym,sym); // TODO: mixed variant