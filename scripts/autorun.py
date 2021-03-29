# Automatically run the relevant chebyshev expansion
# script for a set range of pms coefficients

import generate_new_cheb_data as gen

def main():

	for i in range(51):
		gen.main(0,i,1e-2)
		gen.main(+2,i,1e-2)
		gen.main(-2,i,1e-2)
	
if __name__ == '__main__':
    main()