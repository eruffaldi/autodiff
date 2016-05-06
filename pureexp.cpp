/**
 * Example of Template Expression 
 * See Also: Fast Reverse-Mode Automatic Differentiation using Expression Templates in C++
 *
 * Improvements:
 */
#include <iostream>
#include "math.h"
struct Stack;
// curiously recurryng template pattern
template <class A>
struct Expression {
	using At = A;
	const A & cast() const { return static_cast<const A&>(*this); }
	double value () const { return cast().value(); }
	//void calc_gradient(Stack & cast) const { cast().calc_gradient(stack); }
};

struct adouble: public Expression<adouble> {
	explicit adouble(double d) : val_(d) {}
	adouble() : val_(0) {}
	template <class R>
	adouble & operator= (const Expression<R> & rhs)
	{
		val_ = rhs.value();
		return *this;
	}
	double value() const {return val_; }
private:
	double val_;
};

template <class A>
struct Sin: public Expression< Sin <A> > {
	// ADEPT performs the computation HERE!
	Sin(const Expression<A>& a) : a_(a.cast()), result_(sin(a_.value())) {}
	double value() const { return result_; }
private:
	const A & a_;
	double result_;
};

template <class L,class R>
struct Multiply: public Expression<Multiply<L,R > > {
	Multiply(const Expression<L>& l, const Expression<R>& r) : r_(r.cast()),l_(l.cast()), result_(l_.value()*r_.value()) {}
	double value() const { return result_; }
private:
	const R & r_;
	const L & l_;
	double result_;
};


template <class L,class R>
struct Add: public Expression<Add<L,R > > {
	Add(const Expression<L>& l, const Expression<R>& r) : r_(r.cast()),l_(l.cast()), result_(l_.value()+r_.value()) {}
	double value() const { return result_; }
private:
	const R & r_;
	const L & l_;
	double result_;
};

template <class L,class R>
struct Sub: public Expression<Sub<L,R > > {
	Sub(const Expression<L>& l, const Expression<R>& r) : r_(r.cast()),l_(l.cast()), result_(l_.value()-r_.value()) {}
	double value() const { return result_; }
private:
	const R & r_;
	const L & l_;
	double result_;
};

template <class A> inline
Sin<A> sin(const Expression<A> & a) { return Sin<A>(a); }

template <class A,class B> inline
Multiply<A,B> operator*(const Expression<A>&a,const Expression<B>&b) { return Multiply<A,B>(a,b); }

template <class A,class B> inline
Add<A,B> operator+(const Expression<A>&a,const Expression<B>&b) { return Add<A,B>(a,b); }

template <class A,class B> inline
Sub<A,B> operator-(const Expression<A>&a,const Expression<B>&b) { return Sub<A,B>(a,b); }


int main(int argc, char * argv[])
{
	adouble q(10);
	adouble b(20.0);
	auto r = q*b + sin(b);
	std::cout << r.value() << std::endl;
	return 0;
}
