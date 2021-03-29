# Automatically run the relevant chebyshev expansion
# script for a set range of pms coefficients

import generate_new_cheb_data as gen

def main():

	for i in range(101):
		gen.main(0,i,1e-6)
		gen.main(+2,i,1e-6)
		gen.main(-2,i,1e-6)
	
if __name__ == '__main__':
    main()