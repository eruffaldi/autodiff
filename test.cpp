#include <iostream>
#include "sympp.hpp"

/*
y = 1
mu = 2
sigma = 3
syms y mu sigma real
f = -0.5*((y-mu)/sigma)^2-log(sigma)-0.5*log(2*pi())
df = jacobian(f,{y,mu,sigma});
% [ (2*mu - 2*y)/(2*sigma^2), -(2*mu - 2*y)/(2*sigma^2), (mu - y)^2/sigma^3 - 1/sigma]

% 

fv = subs(f,{y,mu,sigma},{1,2,3});
double(fv) % -2.0731
dfv = subs(df,{y,mu,sigma},{1,2,3});
double(dfv) % 0.1111   -0.1111   -0.2963
*/

void descend(sym v, std::string q = std::string())
{
	std::cout << q << typeid(*v.p_).name() << " nparents " << v.nparents() << " ptr " << v.p_.get() << std::endl;
	q += ' ';
	for(int i = 0; i < v.nparents(); i++)
	{
		descend(v.parent(i),q);
	}
}

sym lognormal(sym y,sym mu,sym sigma)
{
	return -0.5*pow((y-mu)/sigma,2)-log(sigma)-0.5*log(2*pi());
}

int main(int argc, char const *argv[])
{
	// Sten compatibility: grad() adj() 
	sym y = sym("y",1.0);
	sym mu = sym("mu",2.0);
	sym sigma = sym("sigma",3.0);

	sym f = lognormal(y,mu,sigma);
	std::cout << "function " << f << std::endl;

	std::cout << "eval " << f.val() << std::endl;

	std::cout << "descend\n";
	 descend(f);
	/*
	// TODO: evaluate when y,mu,sigma are assigned
	auto fs = subs(f,{y},{var(0.1)});
	std::cout << "function replacement " << fs << std::endl;
	*/
	std::cout << "symbolic jacobian\n";
	jacob J(f,{y,mu,sigma});
	for(int i = 0; i < J.gradients(); i++)
	{
		std::cout << J.gradientVar(i) << ": " << J.gradient(i) << " => " << J.gradient(i).val() << std::endl;
	//	std::cout << "descend\n";
		 //descend(J.gradient(i),"*");
	}
	
	std::cout << "numeric jacobian\n";
	jacobnum J2(f,{y,mu,sigma});
	for(int i = 0; i < J2.gradients(); i++)
	{
		std::cout << J2.gradientVar(i) << ": " << J2.gradient(i) << std::endl;
	}
	return 0;
}