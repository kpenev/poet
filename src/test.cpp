//#include "Common.h"
//#include <iostream>
//#include <solve_polynomial.h>
#include <cstdlib>
#include <iostream>
#include <sstream>

#define STR(x) #x
#define STRING(str) STR(str)
#define TEST_VALUE 1.3

/*void output_array(std::valarray<double> a, std::string name)
{
	std::cout << name << std::endl;
	for(int i=0; i<a.size(); i++)
		std::cout << a[i] << ", ";
	std::cout << std::endl;
}*/

int main()
{
/*	std::valarray<double> exact_solutions=solve_cubic(-0.04, 0.96, 1e-12, 0.0),
		approx_solutions=solve_cubic(-0.04, 0.96, 1, 1e-12);
	output_array(exact_solutions, "quadratic");
	output_array(approx_solutions, "cubic");*/
	std::cout << "test_value=" STRING(TEST_VALUE) << std::endl;
	double *p=NULL;
	std::istringstream is("12,3.2,1");
	while(!is.eof() && is) {
		double v;
		is >> v;
		std::cout << "extracted: " << v << std::endl;
	}
	return 0;
}
