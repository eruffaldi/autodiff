/*
Automatic Differentiation with Matrix Support. 
Emanuele Ruffaldi 2015 end of year project
 */
#include "sympp.hpp"
#include <sstream>

/**
 * Special Scalar Constants with pre-allocated static instances. This reduces the allocation for very common symbols, and, 
 * at the same time allows for nice printing of transcendental numbers (e.g. pi and e).
 *
 */
struct vconstspecial : public isym
{
	/// enumeration of the possibile constants
	enum Const { ZERO, ONE, NEGONE, TWO, HALF, PI, E, LOG2};

	/// constructor with all the inspectionable information
	vconstspecial(Const aid, double avalue, std::string aname): id(aid),value(avalue),name(aname)
	{}

public:
	double value;
	std::string name;
	Const id;

	 std::string sig() const override
	 {
	 	std::ostringstream os; 
	 	os << name;
	 	return os.str();

	 }
	 void print(std::ostream & os) const override { os << name; }	 
	 double operator()() const override { return value; }
	 bool isconst() const override { return true; }
	 std::shared_ptr<isym> diff(int p) const override { return zero.p_; }
	 int nparents() const override { return 0; }
	 std::shared_ptr<isym> parent(int p) const override { return nullptr; }

	 static sym e;
	 static sym pi;
	 static sym log2;
	 static sym two;
	 static sym half;
	 static sym one;
	 static sym negone;
	 static sym zero;
};

sym vconstspecial::e(std::make_shared<vconstspecial>(vconstspecial::E,2.7182818284590452353602874,"e"));
sym vconstspecial::pi(std::make_shared<vconstspecial>(vconstspecial::PI,3.141592653589793238462643,"pi"));
sym vconstspecial::log2(std::make_shared<vconstspecial>(vconstspecial::LOG2,0.693147180559945,"log2"));
sym vconstspecial::two(std::make_shared<vconstspecial>(vconstspecial::TWO,2,"2.0"));
sym vconstspecial::half(std::make_shared<vconstspecial>(vconstspecial::HALF,0.5,"0.5"));
sym vconstspecial::one(std::make_shared<vconstspecial>(vconstspecial::ONE,1,"1.0"));
sym vconstspecial::negone(std::make_shared<vconstspecial>(vconstspecial::NEGONE,-1,"-1.0"));
sym vconstspecial::zero(std::make_shared<vconstspecial>(vconstspecial::ZERO,0,"0.0"));

/**
 * Constant value not belonging to the ones provided by the user
 */
struct vconst : public isym
{
	vconst(double d) : value(d) {}

	double value;

	 std::string sig() const override
	 {
	 	std::ostringstream os; 
	 	os << value;
	 	return os.str();

	 }

	 bool isconst() const override { return true; }
	 void print(std::ostream & os) const override { os << value; }
	 double operator()() const override { return value; }

	 std::shared_ptr<isym> diff(int p) const override { return vconstspecial::zero.p_; }
	 int nparents() const override { return 0; }
	 std::shared_ptr<isym> parent(int p) const override { return nullptr; }
};



/**
 * Variable Symbol with name and specification
 * THIS is partially immutable in the sense that name and spec cannot be changes because expression construction relies on them. 
 * The value can be modified at will for supporting expression evaluation
 */
struct vsymbol : public isym
{
	std::string name;
	double value;

	/// name and scalar value
	vsymbol(std::string n, double v = 0) : name(n) { value = v; }

	/// name and matrix value
	//vsymbol(std::string n, Eigen::MatrixXd v) : name(n),value(v) { spec = matrixspec(v.rows(),v.cols()); }


	 std::string sig() const  override { return name; }

	 void print(std::ostream & os) const override { os << name; }
	 double operator()() const override { return value; }

	 // TODO matrix version

	 // TODO matrix version
	 std::shared_ptr<isym> diff(int p) const override { return vconstspecial::zero.p_; }
	 int nparents() const override { return 0; }
	 std::shared_ptr<isym> parent(int p) const override { return nullptr; }

};

/**
 * Base for unary operations
 *
 * TODO: inherit the spec
 */
struct vunop : public isym
{
	std::shared_ptr<isym> up;

    vunop(std::shared_ptr<isym> a) : up(a) {}
    vunop(sym a) : up(a.p_) {}

    virtual std::string fxop() const = 0;
    std::string sig() const  override { return fxop() + std::string("::") + up->sig(); }
    void print(std::ostream & os) const override { os << fxop() << '('; up->print(os); os << ')'; }

	 int nparents() const override { return 1; }
	 std::shared_ptr<isym> parent(int p) const override { return p == 0 ? up : nullptr; }

};

/**
 * Base for binary operations
 */
struct vbinop : public isym
{	
	std::shared_ptr<isym> first,second;
    
    vbinop(std::shared_ptr<isym> a, std::shared_ptr<isym> b) : first(a),second(b) {}
    vbinop(sym a, sym b) : first(a.p_),second(b.p_) {}
    
    virtual char op() const = 0;
    std::string sig() const  override { return op() + std::string("::") + first->sig() + "::" + second->sig(); }


	 int nparents() const override { return 2; }
	 std::shared_ptr<isym> parent(int p) const override { return p == 0 ? first : p == 1 ? second : nullptr; }
};


/**
 * Summation

 TODO: verification of operands in matrix form, and computation of output matrix form
 */
struct vsumop: public vbinop
{
    using vbinop::vbinop;

    char op() const override { return '+'; }
    void print(std::ostream & os) const override { os << "(" ; first->print(os); os << ")+("; second->print(os); os << ")"; }
    double operator()() const override { return (*first)()+(*second)(); }

	std::shared_ptr<isym> diff(int p) const override { return vconstspecial::one.p_; }
};

/**
 * Difference

 TODO: verification of operands in matrix form, and computation of output matrix form
 */
struct vdiffop: public vbinop
{
    using vbinop::vbinop;
    
    char op() const override { return '-'; }
    void print(std::ostream & os) const override { os << "(" ; first->print(os); os << ")-("; second->print(os); os << ")"; }
    double operator()() const override { return (*first)()-(*second)(); }

	std::shared_ptr<isym> diff(int p) const override { return p == 0 ? vconstspecial::one.p_: vconstspecial::negone.p_; }
};

/**
 * Multiplication 
 *
 TODO: matrix form computation and differentiation
 */
struct vmulop: public vbinop
{
    using vbinop::vbinop;
    
    char op() const override { return '*'; }
    void print(std::ostream & os) const override { os << "(" ; first->print(os); os << ")*("; second->print(os); os << ")"; }
    double operator()() const override { return (*first)()*(*second)(); }

	std::shared_ptr<isym> diff(int p) const override { return p == 0 ? second: first; }
};

/**
 * Division a/b 
 *
 * TODO: matrix form
 */
struct vdivop: public vbinop
{
    using vbinop::vbinop;
    
    char op() const override { return '/'; }
    void print(std::ostream & os) const override { os << "(" ; first->print(os); os << ")/("; second->print(os); os << ")"; }
    double operator()() const override { return (*first)()/(*second)(); }

    // f = x/y
    // df/dx = 1/y == f/x
    // df/dy = -x/y^2 == -f/y
	std::shared_ptr<isym> diff(int p) const override { return (p == 0 ? invert(second): (-first->tosym()*pow(second->tosym(),-2))).p_; }

};

struct vsinop: public vunop
{
	using vunop::vunop;
	std::string fxop() const  override { return "sin"; }
	double operator()() const override { return sin((*up)()); }
	std::shared_ptr<isym> diff(int p) const override { return cos(up->tosym()).p_; }
};

struct vcosop: public vunop
{
    using vunop::vunop;
    std::string fxop() const  override { return "cos"; }
    double operator()() const override { return cos((*up)()); }
 	std::shared_ptr<isym> diff(int p) const override { return (-sin(up->tosym())).p_; }
};

/// TODO: derivative for matrix
struct vexpop: public vunop
{
    using vunop::vunop;
    std::string fxop() const  override { return "exp"; }
    double operator()() const override { return exp((*up)()); }
 	std::shared_ptr<isym> diff(int p) const override { return (exp(up->tosym())).p_; }
};

/// TODO: derivative for matrix
struct vlogop: public vunop
{
    using vunop::vunop;
    std::string fxop() const  override { return "log"; }
    double operator()() const override { return log((*up)()); }
 	std::shared_ptr<isym> diff(int p) const override { return (invert(up->tosym())).p_; }
};

/// TODO: matrix form
struct vpowconstop: public vunop
{
	int power;

	vpowconstop(sym v, int q) : vunop(v), power(q) {}
	vpowconstop(std::shared_ptr<isym> v, int q) : vunop(v), power(q) {}
    std::string fxop() const  override { return "pow"; }
	 std::string sig() const  override { return "pow::" + std::to_string(power) + "_" + up->sig(); }
    double operator()() const override { return pow((*up)(),power); }
    void print(std::ostream & os) const override { os << "pow("; up->print(os); os << "," << power << ")"; }

    // f = x^k
    // df/dx = k x^(k-1)
    //   k=-1 
    //	 k=0 0
    //	 k=1 1
    //	 k=2 2x
    //	 otherwise generic
 	std::shared_ptr<isym> diff(int p) const override {
 		switch(power)
 		{
 			case 0: return zero().p_;
 			case 1: return one().p_;
 			case 2: return (two()*up->tosym()).p_; 
 			default:
 				return (power*pow(up->tosym(),power-1)).p_;
 		}
 	}
};

/// TODO: matrix form
struct vpowop: public vbinop
{
    using vbinop::vbinop;
    char op() const override { return 'a'; }
	 std::string sig() const  override { return "pow::" + first->sig() + "_" + second->sig(); }
    void print(std::ostream & os) const override { os << "pow("; first->print(os); os << ","; second->print(os); os << ")"; }
     double operator()() const override { return pow((*first)(),(*second)()); }

     // f = x^y
     // df/dx = x^(y-1) y = y f / x
     // df/dy = x^y log(x) = f * log(x)
	 std::shared_ptr<isym>  diff(int p) const override
	 {
	 	if(p == 0)
		 	return (second->tosym()*pow(first->tosym(),second->tosym()-one())).p_;
		else
		 	return (tosym() *log(first->tosym())).p_;
	 }

};

sym isym::tosym() const  { return sym(const_cast<isym*>(this)->shared_from_this()); }

sym operator + (const sym & a, const sym & b) {	return sym(std::make_shared<vsumop>(a,b)); }
sym operator - (const sym & a, const sym & b) {	return sym(std::make_shared<vdiffop>(a,b));}
sym operator / (const sym & a, const sym & b) {	return sym(std::make_shared<vdivop>(a,b));}
sym operator * (const sym & a, const sym & b) {	return a == one() ? b : sym(std::make_shared<vmulop>(a,b));}

sym operator + (double a, const sym & b) {	return sym(std::make_shared<vsumop>(b,sym(std::make_shared<vconst>(a)))); }
sym operator - (double a, const sym & b) {	return sym(std::make_shared<vdiffop>(b,sym(std::make_shared<vconst>(a)))); }
sym operator / (double a, const sym & b) {	return sym(std::make_shared<vdivop>(b,sym(std::make_shared<vconst>(a)))); }
sym operator * (double a, const sym & b) {	return sym(std::make_shared<vmulop>(b,sym(std::make_shared<vconst>(a)))); }

sym operator + (const sym & a, double b) {	return sym(std::make_shared<vsumop>(a,sym(std::make_shared<vconst>(b))));	 }
sym operator - (const sym & a, double b) { 	return sym(std::make_shared<vdiffop>(a,sym(std::make_shared<vconst>(b))));  }
sym operator / (const sym & a, double b) {	return sym(std::make_shared<vdivop>(a,sym(std::make_shared<vconst>(b))));	 }
sym operator * (const sym & a, double b) {	return sym(std::make_shared<vmulop>(a,sym(std::make_shared<vconst>(b))));	 }

sym& sym::operator += (const sym & a)  { sym r =  *this + a; p_.swap(r.p_); return *this; }
sym& sym::operator *= (const sym & a)  { sym r =  *this * a; p_.swap(r.p_); return *this; }
sym& sym::operator /= (const sym & a)  { sym r =  *this / a; p_.swap(r.p_); return *this; }
sym& sym::operator -= (const sym & a)  { sym r =  *this - a; p_.swap(r.p_); return *this; }

/// change of sign of a divop flips terms
sym sym::operator - () const  
{
	vdiffop * p = dynamic_cast<vdiffop*>(p_.get());
	if(!p)
		return negone()* *this;
	else
		return sym(std::make_shared<vdiffop>(parent(1),parent(0)));
}

sym two() { return vconstspecial::two; }
sym zero() { return vconstspecial::zero; }
sym one() { return vconstspecial::one; }
sym negone() { return vconstspecial::negone; }
sym e() { return vconstspecial::e; }
sym log2() { return vconstspecial::log2; }
sym half() { return vconstspecial::half; }
sym pi() { return vconstspecial::pi; }

sym log(sym x) { return sym(std::make_shared<vlogop>(x)); }
sym exp(sym x) { return sym(std::make_shared<vexpop>(x)); }
sym sin(sym x) { return sym(std::make_shared<vsinop>(x)); }
sym cos(sym x) { return sym(std::make_shared<vcosop>(x)); }
sym pow(sym x, double q) { return sym(std::make_shared<vpowconstop>(x,q)); }
sym pow(sym x1, sym x2) { return sym(std::make_shared<vpowop>(x1,x2)); }
sym invert(sym x) { return sym(std::make_shared<vpowconstop>(x,-1)); }

// TODO: sqrt of matrices has special semantics
sym sqrt(sym x) { return sym(std::make_shared<vpowconstop>(x,0.5)); }

sym::sym(std::string a, double v): p_(std::make_shared<vsymbol>(a,v))
{
}

/// constant
sym::sym(double d): p_(std::make_shared<vconst>(d))
{

}

/// constant
sym::sym(int d): p_(std::make_shared<vconst>(d))
{

}



#if 0
/// compactifies branches of the expression graph in which all terms are constant
/// for a given non leaf: op(x1,...,xn):
///		all xi are constant => replace with val()
///		for any xi not previously constant => replace it with constant
///
/// REMEMBER sym is immutable so we need to make a clone with replacements
sym sym::solveconst()
{
	#if 0
	class constsolver
	{
	public:
		sym solve(sym x)
		{
			sym r(x);
			double v = 0;
			if(descend(r,v))
				return sym(v); // full const
			else
				return r;
		}



		bool descend(sym & x, double & v)
		{
			int n = x.nparents();
			if(n == 0)
			{
				isym * p = x.p_get();
				if((vconstspecial * ps = dynamic_cast<vconstspecial>(p)) != 0)
				{	
					v = ps->value;
					return true;
				}
				else if((vconst * pv = dynamic_cast<vconst>(p)) != 0)
				{
					v = ps->value;
					return true;
				}	
				else
					return false;
			}
			else
			{
				for(int i = 0; i < n; i++)
				{
					sym xx(x.parent(i));
					double w;
					if(descend(xx,w))
					{
						// need replace one OR all with w
					}
					else
					{

					}
				}
			}
		}
	};


	constsolver cs;
	return cs.solve(*this);	
	#endif
	return *this;
}
#endif

jacob::jacob(sym root,std::initializer_list<sym> syms)
{
	for(auto x: syms)
		targets_.push_back(x.p_);
	adjoint_.emplace(std::pair<std::shared_ptr<isym>, sym >(root.p_, one()));
	descend(root,true);
}

sym jacob::gradient(sym s) const 
{
	auto it = adjoint_.find(s.p_);
	if(it == adjoint_.end())
		return zero();
	else
		return sym(it->second);
}


void jacob::descend(sym v, bool isroot)
{
	int n = v.nparents();

	if(n == 0) //leaf
	{
	}
	else
	{
		auto itm = adjoint_.find(v.p_);
		if(itm == adjoint_.end())
		{
			std::cerr << "error wrong order visit\n";
			return;
		}
		for(int i = 0; i < n; i++)
		{
			sym p = v.parent(i);
			auto it = adjoint_.find(p.p_);

			/// TODO: matrix form has different way of doing this
			sym r = (v.diff(i) * itm->second);
			if(it == adjoint_.end())
				adjoint_.emplace(std::pair<std::shared_ptr<isym>, sym >(p.p_, r));
			else
				it->second += r;
			descend(p,false);
		}
	}
}




jacobnum::jacobnum(sym root,std::initializer_list<sym> syms): root_(root)
{
	for(auto x: syms)
		targets_.push_back(x.p_);
	descendG(root_);
	update();
}

void jacobnum::update()
{
	adjoint_.clear();
	adjoint_.emplace(std::pair<std::shared_ptr<isym>, double >(root_.p_, 1));
	descendR(root_);
}

double jacobnum::gradient(sym s) const 
{
	auto it = adjoint_.find(s.p_);
	if(it == adjoint_.end())
		return 0;
	else
		return it->second;
}


void jacobnum::descendR(sym v)
{
	int n = v.nparents();

	if(n == 0) //leaf
	{
	}
	else
	{
		auto itm = adjoint_.find(v.p_);
		if(itm == adjoint_.end())
		{
			std::cerr << "error wrong order visit\n";
			return;
		}
		auto itG = diff_.find(v.p_);
		for(int i = 0; i < n; i++)
		{
			double r = itG->second[i].val()*itm->second;
			sym p = v.parent(i);
			auto it = adjoint_.find(p.p_);
			if(it == adjoint_.end())
				adjoint_.emplace(std::pair<std::shared_ptr<isym>, double >(p.p_, r));
			else
				it->second += r;
			descendR(p);
		}
	}
}

void jacobnum::descendG(sym v)
{
	int n = v.nparents();
	if(n == 0) //leaf
	{
	}
	else
	{
		{
		auto itG = diff_.find(v.p_);
		if(itG != diff_.end())
			return;
		}
		{
			std::vector<sym> w;
			w.reserve(n);
			for(int i = 0; i < n; i++)
				w.push_back(v.diff(i));
			diff_[v.p_].swap(w);
		}
		for(int i = 0; i < n; i++)
			descendG(v.parent(i));
	}
}
