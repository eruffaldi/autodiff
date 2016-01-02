/**
 http://www2.maths.ox.ac.uk/~gilesm/files/NA-08-01.pdf

 Matrix Operations:

 Addition: C=A+B
 	aA=aC
 	aB=aC
 Product:  C=AB
 	aA=aC B'
 	aB=A' aC
 Inverse:  C=inv A
 	aA=-C' aC C'
 Det: 
 	aA= aC C inv(A')
 Induced Special: C = (inv A) B
 	aB = inv(A') aC
 	aA = -aB C'
 Induced Special: C = B'AB
 	aA = B aC B'
 	aB = A B aC' + A'B aC
 Induced Special: C = B'invA B
 	aA = - invA' B aC B' invA'
 	aB = invA B aC' + invA' B aC

 Frobenius Norm:
 	aA = aB inv(B) A
 Spectral/Euclidean Norm:
 	aA = aB U1 V1'   with first columns of SVD 
 Trace? 	 
 	Maybe: http://www.tc.umn.edu/~nydic001/docs/unpubs/Schonemann_Trace_Derivatives_Presentation.pdf
			http://web.stanford.edu/~jduchi/projects/matrix_prop.pdf
 	C = Tr(A) = Sum Aii
 	BUT
	 	C = Sum lambda_i 
	 	det A = Prod lambda_i
 	SO
 		the trace is the derivative of the determinant
	
	1) dC = d(f/A) dA = A dA
		note that for C=detA we had  dC = detA Tr(inv A dA)
	2) Tr(aC' dC) = Tr(aA' dA) 
	   Tr(aC' A dA) = Tr(aA' dA)
	   so
 	 	   	aC' A = aA'  =>   aA = A' aC

 	REALLY?



 Also: polynomial and exponential

 */
#include "msympp.hpp"



msym imsym::tosym() const  { return msym(const_cast<imsym*>(this)->shared_from_this()); }



#if 0
struct vconstmat : public imsym
{
	//virtual std::shared_ptr<isym> diff(int index) const = 0;
	 std::string sig() const override { return ""; }
	 void print(std::ostream & os) const override
	{

	}

	 Eigen::MatrixXd operator()() const override
	 {

	 }

	 std::shared_ptr<imsym> diff(int p) const override
	 {

	 }

	 int nparents() const override { return 0; }
	 std::shared_ptr<imsym> parent(int p) const override  {}
	 msym parentadjointS(int i) const override {}
	 Eigen::MatrixXd  parentadjointN(int i) const override {}
};
#endif

/**
 * Constant Matrix
 */
struct mvconstmat : public imsym
{
	mvconstmat(Eigen::MatrixXd value): value_(value),spec_(value.rows(),value.cols()) {}

	//virtual std::shared_ptr<isym> diff(int index) const = 0;
	std::string sig() const override { return "vconstmat"; }

	void print(std::ostream & os) const override
	{
		os << "vconstmat " << value_;
	}

	 Eigen::MatrixXd operator()() const override
	 {
	 	return value_;
	 }

	 std::shared_ptr<imsym> diff(int p) const override
	 {
	 	return nullptr; // TODO: zero with size
	 }

		bool isconst() const override { return true; }
	 int nparents() const override { return 0; }
	 std::shared_ptr<imsym> parent(int p) const override { return nullptr; }
	 msym parentadjointS(int i) const override { return nullptr; }
	 Eigen::MatrixXd  parentadjointN(int i) const override { return Eigen::MatrixXd::Zero(rows(),cols()); }

	 Eigen::MatrixXd value_;
};

/*
 * Matrix array of symbols (valuable)
 */
struct mvsymbol : public imsym
{
	mvsymbol(std::string name, matrixspec spec) {}

	//virtual std::shared_ptr<isym> diff(int index) const = 0;
	 std::string sig() const override { return "mvsymbol"; }

	void print(std::ostream & os) const override
	{
		os << "mvsymbol ";
	}

	 Eigen::MatrixXd operator()() const override
	 {
	 	return value_;
	 }

	 std::shared_ptr<imsym> diff(int p) const override
	 {
	 	return nullptr; // TODO
	 }

	 bool isconst() const override { return false; }
	 int nparents() const override { return 0; }
	 std::shared_ptr<imsym> parent(int p) const override { return nullptr; }
	 msym parentadjointS(int i) const override { return nullptr; }
	 Eigen::MatrixXd  parentadjointN(int i) const override { return Eigen::MatrixXd::Zero(rows(),cols()); }
	 bool set(const Eigen::MatrixXd &x) override { value_ = x;  return true; }

	 Eigen::MatrixXd value_;
};

/**
 * Base for binary operations
 */
struct mvbinop : public isym
{	
	std::shared_ptr<imsym> first,second;
    
    mvbinop(std::shared_ptr<imsym> a, std::shared_ptr<imsym> b) : first(a),second(b) {}
    mvbinop(msym a, msym b) : first(a.p_),second(b.p_) {}
    
    virtual char op() const = 0;
    std::string sig() const  override { return op() + std::string("::") + first->sig() + "::" + second->sig(); }

	 int nparents() const override { return 2; }
	 std::shared_ptr<imsym> parent(int p) const override { return p == 0 ? first : p == 1 ? second : nullptr; }
};

struct mvunop : public imsym
{
	std::shared_ptr<imsym> up;

    vunop(std::shared_ptr<imsym> a) : up(a) {}
    vunop(msym a) : up(a.p_) {}

    virtual std::string fxop() const = 0;
    std::string sig() const  override { return fxop() + std::string("::") + up->sig(); }
    void print(std::ostream & os) const override { os << fxop() << '('; up->print(os); os << ')'; }

	 int nparents() const override { return 1; }
	 std::shared_ptr<imsym> parent(int p) const override { return p == 0 ? up : nullptr; }

};

struct mvsumop: public mvbinop
{
    using mvbinop::mvbinop;

    char op() const override { return '+'; }
    void print(std::ostream & os) const override { os << "(" ; first->print(os); os << ")+("; second->print(os); os << ")"; }
    Eigen::MatrixXd operator()() const override { return (*first)()+(*second)(); }

	//std::shared_ptr<isym> diff(int p) const override { return vconstspecial::one.p_; }
};

struct mvdiffop: public mvbinop
{
    using mvbinop::mvbinop;

    char op() const override { return '-'; }
    void print(std::ostream & os) const override { os << "(" ; first->print(os); os << ")-("; second->print(os); os << ")"; }
    Eigen::MatrixXd operator()() const override { return (*first)()-(*second)(); }

	//std::shared_ptr<isym> diff(int p) const override { return vconstspecial::one.p_; }
};

struct mvmulop: public mvbinop
{
    using mvbinop::mvbinop;
    
    char op() const override { return '*'; }
    void print(std::ostream & os) const override { os << "(" ; first->print(os); os << ")*("; second->print(os); os << ")"; }
    Eigen::MatrixXd operator()() const override { return (*first)()*(*second)(); }

	//std::shared_ptr<isym> diff(int p) const override { return p == 0 ? second: first; }
};

struct mvinvert: public mvunop
{
    using mvunop::mvunop;
    std::string fxop() const  override { return "invert"; }
    Eigen::MatrixXd operator()() const override { return (*up)().inverse(); }
 	
 	//std::shared_ptr<isym> diff(int p) const override { return (invert(up->tosym())).p_; }
};



#if 0

/// TODO: evaluation of transpose
/// TODO: evaluation of print for displaying the transposed form
struct vtranspose: public vunop
{
    using vunop::vunop;
	vtranspose(std::shared_ptr<isym> v) : vunop(v) { spec = v->spec.transpose(); }

    std::string fxop() const  override { return "transpose"; }
    double operator()() const override { return NAN; }

    // TODO NOT IMPLEMENTED
 	std::shared_ptr<isym> diff(int p) const override { return tosym().p_; }
};

/// TODO: evaluation and diff of det
struct vdet: public vunop
{
    using vunop::vunop;
	vdet(std::shared_ptr<isym> v) : vunop(v) { spec = matrixspec(); }

    std::string fxop() const  override { return "det"; }
    double operator()() const override { return NAN; }

    // TODO NOT IMPLEMENTED
 	std::shared_ptr<isym> diff(int p) const override { return tosym().p_; }
};

/// TODO: evaluation and diff of trace
struct vtrace: public vunop
{
    using vunop::vunop;
	vtrace(std::shared_ptr<isym> v) : vunop(v) { spec = matrixspec(); }

    std::string fxop() const  override { return "trace"; }
    double operator()() const override { return NAN; }

    // TODO NOT IMPLEMENTED
 	std::shared_ptr<isym> diff(int p) const override { return tosym().p_; }
};

/// TODO: evaluation and diff of frobenius
struct vfrobenius: public vunop
{
    using vunop::vunop;
	vfrobenius(std::shared_ptr<isym> v) : vunop(v) { spec = matrixspec(); }

    std::string fxop() const  override { return "frobenius"; }

    /// sqrt(trace(A A'))
    double operator()() const override { return NAN; }

    // TODO NOT IMPLEMENTED
 	std::shared_ptr<isym> diff(int p) const override { return tosym().p_; }
};


sym::sym(std::string name, const Eigen::MatrixXd & value) : p_(std::make_shared<vsymbol>(name,value))
{

}
sym transpose(sym x) { if(x.spec().isscalar()) return x; else return sym(std::make_shared<vtranspose>(x)); }
sym det(sym x) { if(x.spec().isscalar()) return x; else return sym(std::make_shared<vdet>(x)); }
sym trace(sym x) { if(x.spec().isscalar()) return x; else return sym(std::make_shared<vtrace>(x)); }
sym frobenius(sym x) { if(x.spec().isscalar()) return x; else return sym(std::make_shared<vfrobenius>(x)); }

#endif


/// variable
msym::msym(std::string name, int rows, int cols, MatrixStorage s = MatrixStorage::Scalar,MatrixContent cc = MatrixContent::Generic ):
	msym(name,matrixspec(rows,cols,s,cc))
{
	
}

/// variable
msym::msym(std::string name, matrixspec spec) : p_(std::make_shared<mvsymbol>(name,spec))
{
	
}


/// valued matrix
msym::msym(const Eigen::MatrixXd & m) : p_(std::make_shared<mvconst>(m))
{	
}


