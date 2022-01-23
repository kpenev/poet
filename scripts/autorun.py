# Automatically run the relevant chebyshev expansion
# script for a set range of pms coefficients

import generate_new_cheb_data as gen

def main():

	for i in range(10):#101):
		gen.main(0,i+0,2e-9)
		gen.main(+2,i+0,2e-9)
		gen.main(-2,i+0,2e-9)
	
if __name__ == '__main__':
    main()