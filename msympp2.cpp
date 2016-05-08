/*
Automatic Differentiation with Matrix Support. 
Emanuele Ruffaldi 2015 end of year project
 */
#include "msympp2.hpp"
#include <sstream>

#if 0
/**
 * Special Scalar Constants with pre-allocated static instances. This reduces the allocation for very common msymbols, and, 
 * at the same time allows for nice printing of transcendental numbers (e.g. pi and e).
 *
 */
struct vconstspecial : public imsym
{
	/// enumeration of the possibile constants
	enum Const { ZERO, ONE, NEGONE, TWO, HALF, PI, E, LOG2};

	/// constructor with all the inspectionable information
	vconstspecial(Const aid, double avalue, std::string aname): id(aid),value(avalue),name(aname)
	{}
	
	bool set(mat_t v)  override { return false; }


public:
	mat_t value;
	std::string name;
	Const id;

	 std::string sig() const override
	 {
	 	std::ostringstream os; 
	 	os << name;
	 	return os.str();

	 }
	 void print(std::ostream & os) const override { os << name; }	 
	 mat_t operator()() const override { return value; }
	 bool isconst() const override { return true; }
	 std::shared_ptr<imsym> diff(int p) const override { return zero.p_; }
	 int nparents() const override { return 0; }
	 std::shared_ptr<imsym> parent(int p) const override { return nullptr; }

	 static msym e;
	 static msym pi;
	 static msym log2;
	 static msym two;
	 static msym half;
	 static msym one;
	 static msym negone;
	 static msym zero;
};

msym vconstspecial::e(std::make_shared<vconstspecial>(vconstspecial::E,2.7182818284590452353602874,"e"));
msym vconstspecial::pi(std::make_shared<vconstspecial>(vconstspecial::PI,M_PI,"pi"));
msym vconstspecial::log2(std::make_shared<vconstspecial>(vconstspecial::LOG2,0.693147180559945,"log2"));
msym vconstspecial::two(std::make_shared<vconstspecial>(vconstspecial::TWO,2,"2.0"));
msym vconstspecial::half(std::make_shared<vconstspecial>(vconstspecial::HALF,0.5,"0.5"));
msym vconstspecial::one(std::make_shared<vconstspecial>(vconstspecial::ONE,1,"1.0"));
msym vconstspecial::negone(std::make_shared<vconstspecial>(vconstspecial::NEGONE,-1,"-1.0"));
msym vconstspecial::zero(std::make_shared<vconstspecial>(vconstspecial::ZERO,0,"0.0"));

// TODO: identity and zero
#endif

/**
 * Constant value not belonging to the ones provided by the user
 */
struct vconst : public imsym
{
	vconst(mat_t d) : value(d) {}

	vconst(double d) : value(d) {}

	mat_t value;

	 std::string sig() const override
	 {
	 	std::ostringstream os; 
	 	os << value;
	 	return os.str();

	 }

	bool set(mat_t v) override { return false; }

	 bool isconst() const override { return true; }
	 void print(std::ostream & os) const override { os << value; }
	 mat_t operator()() const override { return value; }

	 std::shared_ptr<imsym> diff(int p) const override { return vconstspecial::zero.p_; }
	 int nparents() const override { return 0; }
	 std::shared_ptr<imsym> parent(int p) const override { return nullptr; }
};




/**
 * Variable msymbol with name and specification
 * THIS is partially immutable in the sense that name and spec cannot be changes because expression construction relies on them. 
 * The value can be modified at will for supporting expression evaluation
 */
struct vmsymbol : public imsym
{
	std::string name;
	mat_t value;

	/// name and scalar value
	vmsymbol(std::string n, mat_t v) : name(n) { value = v; }

	vmsymbol(std::string n, double v = 0) : name(n) { value = v; }

	bool set(mat_t v) override { value = v; return true; }

	bool set(double v) override { value = v; return true; }

	/// name and matrix value
	//vmsymbol(std::string n, Eigen::MatrixXd v) : name(n),value(v) { spec = matrixspec(v.rows(),v.cols()); }


	 std::string sig() const  override { return name; }

	 void print(std::ostream & os) const override { os << name; }
	 mat_t operator()() const override { return value; }

	 // TODO matrix version

	 // TODO matrix version
	 std::shared_ptr<imsym> diff(int p) const override { return vconstspecial::zero.p_; }
	 int nparents() const override { return 0; }
	 std::shared_ptr<imsym> parent(int p) const override { return nullptr; }

};

#if 0
/**
 * Base for unary operations
 *
 * TODO: inherit the spec
 */
struct vunop : public imsym
{
	std::shared_ptr<imsym> up;

    vunop(std::shared_ptr<imsym> a) : up(a) {}
    vunop(msym a) : up(a.p_) {}

    virtual std::string fxop() const = 0;
	bool set(mat_t v) override { return false; }
    std::string sig() const  override { return fxop() + std::string("::") + up->sig(); }
    void print(std::ostream & os) const override { os << fxop() << '('; up->print(os); os << ')'; }

	 int nparents() const override { return 1; }
	 std::shared_ptr<imsym> parent(int p) const override { return p == 0 ? up : nullptr; }

};

/**
 * Base for unary operations
 *
 * TODO: inherit the spec
 */
struct vunscalarop : public imsym
{
	std::shared_ptr<imsym> up;

    vunscalarop(std::shared_ptr<imsym> a) : up(a) {}
    vunscalarop(msym a) : up(a.p_) {}

    virtual std::string fxop() const = 0;
	bool set(mat_t v) override { return false; }
    std::string sig() const  override { return fxop() + std::string("::") + up->sig(); }
    void print(std::ostream & os) const override { os << fxop() << '('; up->print(os); os << ')'; }

	 int nparents() const override { return 1; }
	 std::shared_ptr<imsym> parent(int p) const override { return p == 0 ? up : nullptr; }

};

/**
 * Base for binary operations
 */
struct vbinop : public imsym
{	
	std::shared_ptr<imsym> first,second;
    
    vbinop(std::shared_ptr<imsym> a, std::shared_ptr<imsym> b) : first(a),second(b) {}
    vbinop(msym a, msym b) : first(a.p_),second(b.p_) {}
    
    virtual char op() const = 0;
	bool set(mat_t v) override { return false; }
    std::string sig() const  override { return op() + std::string("::") + first->sig() + "::" + second->sig(); }


	 int nparents() const override { return 2; }
	 std::shared_ptr<imsym> parent(int p) const override { return p == 0 ? first : p == 1 ? second : nullptr; }
};

/**
 * Base for binary operations
 */
struct vbinscalarop : public imsym
{	
	std::shared_ptr<imsym> first,second;
    
    vbinscalarop(std::shared_ptr<imsym> a, std::shared_ptr<imsym> b) : first(a),second(b) {}
    vbinscalarop(msym a, msym b) : first(a.p_),second(b.p_) {}
    
    virtual char op() const = 0;
	bool set(mat_t v) override { return false; }
    std::string sig() const  override { return op() + std::string("::") + first->sig() + "::" + second->sig(); }


	 int nparents() const override { return 2; }
	 std::shared_ptr<imsym> parent(int p) const override { return p == 0 ? first : p == 1 ? second : nullptr; }
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
    mat_t operator()() const override { return (*first)()+(*second)(); }

	std::shared_ptr<imsym> diff(int p) const override { return vconstspecial::one.p_; }
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
    mat_t operator()() const override { return (*first)()-(*second)(); }

	std::shared_ptr<imsym> diff(int p) const override { return p == 0 ? vconstspecial::one.p_: vconstspecial::negone.p_; }
};

/**
 * Memberwise Multiplication 
 *
 */
struct vmulscalarop: public vbinscalarop
{
    using vbinscalarop::vbinscalarop;
    
    char op() const override { return '*'; }
    void print(std::ostream & os) const override { os << "(" ; first->print(os); os << ")*("; second->print(os); os << ")"; }
    mat_t operator()() const override { return (*first)().array()*(*second)().array(); }

    /// out[i,j] = a[i,j] * b[i,j]
	mat_m parentadjointN(int i, const mat_m & adjoint) const override
	{ 
		return i == 0 ? (*first)() : (*second)();
	}
};

/**
 * Multiplication: G(X)*H(X)
 */
struct vmulop: public vbinop
{
    using vbinop::vbinop;
    
    char op() const override { return '*'; }
    void print(std::ostream & os) const override { os << "(" ; first->print(os); os << ")*("; second->print(os); os << ")"; }
    mat_t operator()() const override { return (*first)()*(*second)(); }

 	// Y = G[k,w] H[w,l]
 	// (I_k kron G'(X)) d(H(X),X) + (H(X) kron I_l) d(G(X),X)
	mat_m parentadjointN(int i, const mat_m & adjoint) const override
	{ 
		if(i == 0)
		{
			auto vfirst = *first();
			return kron(stridematrix(vfirst.rows()),vfirst.transpose);
		}
		else
		{
			auto vsecond = *second();
			return kron(vsecond,stridematrix(vsecond.cols()));
		}
	}
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
    mat_t operator()() const override { return (*first)()/(*second)(); }

    // f = x/y
    // df/dx = 1/y == f/x
    // df/dy = -x/y^2 == -f/y
	std::shared_ptr<imsym> diff(int p) const override { return (p == 0 ? invert(second): (-first->tomsym()*pow(second->tomsym(),-2))).p_; }

};

struct vsinop: public vunscalarop
{
	using vunscalarop::vunscalarop;
	std::string fxop() const  override { return "sin"; }
	mat_t operator()() const override { return sin((*up)()); }
	std::shared_ptr<imsym> diff(int p) const override { return cos(up->tomsym()).p_; }
};

struct vcosop: public vunscalarop
{
    using vunscalarop::vunscalarop;
    std::string fxop() const  override { return "cos"; }
    mat_t operator()() const override { return cos((*up)()); }
 	std::shared_ptr<imsym> diff(int p) const override { return (-sin(up->tomsym())).p_; }
};

struct vlogscalarop: public vunscalarop
{
    using vunscalarop::vunscalarop;
    std::string fxop() const  override { return "log"; }
    mat_t operator()() const override { return log((*up)()); }
 	std::shared_ptr<imsym> diff(int p) const override { return (invert(up->tomsym())).p_; }
};

struct vexpscalarop: public vunscalarop
{
    using vunscalarop::vunscalarop;
    std::string fxop() const  override { return "exp"; }
    mat_t operator()() const override { return exp((*up)()); }
 	std::shared_ptr<imsym> diff(int p) const override { return (exp(up->tomsym())).p_; }
};

struct vexpop: public vunop
{
    using vunop::vunop;
    std::string fxop() const  override { return "exp"; }
    mat_t operator()() const override { return exp((*up)()); }
 	std::shared_ptr<imsym> diff(int p) const override { return (exp(up->tomsym())).p_; }
};


struct vlogop: public vunop
{
    using vunop::vunop;
    std::string fxop() const  override { return "log"; }
    mat_t operator()() const override { return log((*up)()); }
 	std::shared_ptr<imsym> diff(int p) const override { return (invert(up->tomsym())).p_; }
};

struct vpowscalarconstop: public vunscalarop
{
	int power;

	vpowscalarconstop(msym v, int q) : vunop(v), power(q) {}
	vpowscalarconstop(std::shared_ptr<imsym> v, int q) : vunop(v), power(q) {}
    std::string fxop() const  override { return "pow"; }
	 std::string sig() const  override { return "pow::" + std::to_string(power) + "_" + up->sig(); }
    mat_t operator()() const override { return pow((*up)(),power); }
    void print(std::ostream & os) const override { os << "pow("; up->print(os); os << "," << power << ")"; }

    // f = x^k
    // df/dx = k x^(k-1)
    //   k=-1 
    //	 k=0 0
    //	 k=1 1
    //	 k=2 2x
    //	 otherwise generic
 	std::shared_ptr<imsym> diff(int p) const override {
 		switch(power)
 		{
 			case 0: return zero().p_;
 			case 1: return one().p_;
 			case 2: return (two()*up->tomsym()).p_; 
 			default:
 				return (power*pow(up->tomsym(),power-1)).p_;
 		}
 	}
};

/// TODO: matrix form
struct vpowscalarop: public vbinscalarop
{
    using vbinscalarop::vbinscalarop;
    char op() const override { return 'a'; }
	 std::string sig() const  override { return "pow::" + first->sig() + "_" + second->sig(); }
    void print(std::ostream & os) const override { os << "pow("; first->print(os); os << ","; second->print(os); os << ")"; }
     mat_t operator()() const override { return pow((*first)(),(*second)()); }

     // f = x^y
     // df/dx = x^(y-1) y = y f / x
     // df/dy = x^y log(x) = f * log(x)
	 std::shared_ptr<imsym>  diff(int p) const override
	 {
	 	if(p == 0)
		 	return (second->tomsym()*pow(first->tomsym(),second->tomsym()-one())).p_;
		else
		 	return (tomsym() *log(first->tomsym())).p_;
	 }

};

msym imsym::tomsym() const  { return msym(const_cast<imsym*>(this)->shared_from_this()); }

msym operator + (const msym & a, const msym & b) {	return msym(std::make_shared<vsumop>(a,b)); }
msym operator - (const msym & a, const msym & b) {	return msym(std::make_shared<vdiffop>(a,b));}
msym operator / (const msym & a, const msym & b) {	return msym(std::make_shared<vdivop>(a,b));}
msym operator * (const msym & a, const msym & b) {	return a == one() ? b : msym(std::make_shared<vmulop>(a,b));}

msym operator + (double a, const msym & b) {	return msym(std::make_shared<vsumop>(b,msym(std::make_shared<vconst>(a)))); }
msym operator - (double a, const msym & b) {	return msym(std::make_shared<vdiffop>(b,msym(std::make_shared<vconst>(a)))); }
msym operator / (double a, const msym & b) {	return msym(std::make_shared<vdivop>(b,msym(std::make_shared<vconst>(a)))); }
msym operator * (double a, const msym & b) {	return msym(std::make_shared<vmulop>(b,msym(std::make_shared<vconst>(a)))); }

msym operator + (const msym & a, double b) {	return msym(std::make_shared<vsumop>(a,msym(std::make_shared<vconst>(b))));	 }
msym operator - (const msym & a, double b) { 	return msym(std::make_shared<vdiffop>(a,msym(std::make_shared<vconst>(b))));  }
msym operator / (const msym & a, double b) {	return msym(std::make_shared<vdivop>(a,msym(std::make_shared<vconst>(b))));	 }
msym operator * (const msym & a, double b) {	return msym(std::make_shared<vmulop>(a,msym(std::make_shared<vconst>(b))));	 }

msym operator + (msym::mat_t a, const msym & b) {	return msym(std::make_shared<vsumop>(b,msym(std::make_shared<vconst>(a)))); }
msym operator - (msym::mat_t a, const msym & b) {	return msym(std::make_shared<vdiffop>(b,msym(std::make_shared<vconst>(a)))); }
msym operator * (msym::mat_t a, const msym & b) {	return msym(std::make_shared<vmulop>(b,msym(std::make_shared<vconst>(a)))); }

msym operator + (const msym & a, msym::mat_t b) {	return msym(std::make_shared<vsumop>(a,msym(std::make_shared<vconst>(b))));	 }
msym operator - (const msym & a, msym::mat_t b) { 	return msym(std::make_shared<vdiffop>(a,msym(std::make_shared<vconst>(b))));  }
msym operator * (const msym & a, msym::mat_t b) {	return msym(std::make_shared<vmulop>(a,msym(std::make_shared<vconst>(b))));	 }

msym& msym::operator += (const msym & a)  { msym r =  *this + a; p_.swap(r.p_); return *this; }
msym& msym::operator *= (const msym & a)  { msym r =  *this * a; p_.swap(r.p_); return *this; }
msym& msym::operator /= (const msym & a)  { msym r =  *this / a; p_.swap(r.p_); return *this; }
msym& msym::operator -= (const msym & a)  { msym r =  *this - a; p_.swap(r.p_); return *this; }

/// change of sign of a divop flips terms
msym msym::operator - () const  
{
	vdiffop * p = dynamic_cast<vdiffop*>(p_.get());
	if(!p)
		return negone()* *this;
	else
		return msym(std::make_shared<vdiffop>(parent(1),parent(0)));
}

msym two() { return vconstspecial::two; }
msym zero() { return vconstspecial::zero; }
msym one() { return vconstspecial::one; }
msym negone() { return vconstspecial::negone; }
msym msym_e() { return vconstspecial::e; }
msym log2() { return vconstspecial::log2; }
msym half() { return vconstspecial::half; }
msym msym_pi() { return vconstspecial::pi; }

msym logdet(msym x) { return msym(std::make_shared<vlogdet>(x)); }

msym trace(msym x) { return msym(std::make_shared<vtrace>(x)); }

msym log(msym x) { return msym(std::make_shared<vlogop>(x)); }
msym exp(msym x) { return msym(std::make_shared<vexpop>(x)); }
msym sin(msym x) { return msym(std::make_shared<vsinop>(x)); }
msym cos(msym x) { return msym(std::make_shared<vcosop>(x)); }
msym pow(msym x, double q) { return msym(std::make_shared<vpowconstop>(x,q)); }
msym pow(msym x1, msym x2) { return msym(std::make_shared<vpowop>(x1,x2)); }
msym invert(msym x) { return msym(std::make_shared<vinvertop>(x)); }

// TODO: sqrt of matrices has special semantics
msym sqrt(msym x) { return msym(std::make_shared<vpowconstop>(x,0.5)); }
#endif

msym::msym(std::string a, double v): p_(std::make_shared<vmsymbol>(a,v))
{
}

msym::msym(std::string a, mat_t v): p_(std::make_shared<vmsymbol>(a,v))
{
}

/// constant
msym::msym(double d): p_(std::make_shared<vconst>(d))
{

}

/// constant
msym::msym(mat_t d): p_(std::make_shared<vconst>(d))
{

}

/// constant
msym::msym(int d): p_(std::make_shared<vconst>(d))
{

}


jacobmsym::jacobmsym(msym root,std::initializer_list<msym> msyms)
{
	for(auto x: msyms)
		targets_.push_back(x.p_);
	adjoint_.emplace(std::pair<std::shared_ptr<imsym>, msym >(root.p_, one()));
	descend(root,true);
}

msym jacobmsym::gradient(msym s) const 
{
	auto it = adjoint_.find(s.p_);
	if(it == adjoint_.end())
		return zero();
	else
		return msym(it->second);
}


void jacobmsym::descend(msym v, bool isroot)
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
			msym p = v.parent(i);
			auto it = adjoint_.find(p.p_);

			/*msym r = (v.diff(i) * itm->second);
			if(it == adjoint_.end())
				adjoint_.emplace(std::pair<std::shared_ptr<imsym>, msym >(p.p_, r));
			else
				it->second += r;
			descend(p,false);
			*/
		}
	}
}

jacobnum::jacobnum(msym root,std::initializer_list<msym> msyms): root_(root)
{
	for(auto x: msyms)
		targets_.push_back(x.p_);
	descendG(root_);
	update();
}

/**
 * During the second evaluation
 */
void jacobnum::update()
{
	pti it;
	if(adjoint_.empty())
	{
		it = adjoint_.emplace(std::pair<std::shared_ptr<imsym>, mat_t >(root_.p_, 1)).first;
	}
	else
	{
		// cleare previous ones to zero
		for(auto & p: adjoint_)
			p.second = 0;
		// re-initialize top most adjoint to 1
		it = adjoint_.find(root_.p_);
		it->second = 1;
	}
	descendR(root_,it);
}

mat_t jacobnum::gradient(msym s) const 
{
	auto it = adjoint_.find(s.p_);
	if(it == adjoint_.end())
		return 0;
	else
		return it->second;
}

void jacobnum::descendR(msym v,pti itm)
{
	int n = v.nparents();

	if(n == 0) //leaf
	{
		return ;
	}
	else
	{
		auto itG = diff_.find(v.p_);
		for(int i = 0; i < n; i++)
		{
			// product between the derivative and the adjoint
			mat_t r = itG->second[i].val()*itm->second;
			msym p = v.parent(i);
			auto it = adjoint_.find(p.p_);
			if(it == adjoint_.end()) // need to create adjoint
				it = adjoint_.emplace(std::pair<std::shared_ptr<imsym>, mat_t >(p.p_, r)).first;
			else
				it->second += r;
			descendR(p,it);
		}
	}
}

void jacobnum::descendG(msym v)
{
	int n = v.nparents();
	if(n == 0) //leaf
	{
	}
	else
	{
		// expression re-use
		{
			auto itG = diff_.find(v.p_);
			if(itG != diff_.end())
				return;
		}
		{
			std::vector<msym> w;
			w.reserve(n);
			for(int i = 0; i < n; i++)
				w.push_back(v.diff(i));
			diff_[v.p_].swap(w);
		}
		for(int i = 0; i < n; i++)
			descendG(v.parent(i));
	}
}
