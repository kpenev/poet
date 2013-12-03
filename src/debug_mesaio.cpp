#include "MESAIO.h"

int main(int, char **argv)
{
	try {
		std::cout << "Generating MESA evolution." << std::endl;
		MESA::Evolution evolution("MESA/tracks");
		std::cout << "Done with MESA evolution" << std::endl;
		const EvolvingStellarQuantity
			*radius=evolution.interpolate_radius(1.0);
		std::cout.setf(std::ios_base::scientific);
		std::cout.precision(16);
		for(double age=1e-4; age<10; age*=1.01)
			std::cout << std::setw(25) << age << std::setw(25) << 
				(*radius)(age) << std::endl;
		return 0;
	} catch(Error::General &ex) {
		std::cout << ex.what() << ": " << ex.get_message() << std::endl;
		return 1;
	}
}
