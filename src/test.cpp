//#include "Common.h"
//#include <iostream>
//#include <solve_polynomial.h>
#include <cstdlib>
#include <iostream>
#include <sstream>

#define STR(x) #x
#define STRING(str) STR(str)
#define TEST_VALUE 1.3

class Base {
public:
	Base() {}
	virtual double operator[](size_t i)=0;
	virtual ~Base() {};
};

class Derived : public Base {
private:
	double *data;
public:
	Derived() {data=new double[50];}
	double operator[](size_t i) {return data[i];}
	~Derived() {delete[] data;}
};

Base &get_base()
{
	return *(new Derived);
}

int main()
{
	Base &base=get_base();
	delete &base;
	return 0;
}
