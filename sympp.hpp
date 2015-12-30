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
	virtual std::shared_ptr<isym> diff(int p) const = 0;
	virtual int nparents() const = 0;
	virtual std::shared_ptr<isym> parent(int p) const = 0;

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