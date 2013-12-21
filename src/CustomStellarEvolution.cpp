/**\file
 *
 * \brief The definitions of the CustomStellarEvolution::Evolution methods.
 */

#include "CustomStellarEvolution.h"

namespace CustomStellarEvolution {

	Evolution::Evolution(const std::string &filename,
			const std::vector<Columns> &format,
			const std::vector<double> &smoothing,
			const std::vector<int> &nodes)
	{
		std::ifstream track(filename.c_str(), std::ios_base::in);
		std::vector<int> column_numbers(NUM_TRACK_QUANTITIES, -1);
		for(size_t i=0; i<NUM_TRACK_QUANTITIES; ++i)
			column_numbers[format[i]]=i;
		std::valarray< std::list<double> >
			quantities=parse_columns(track, column_numbers);
		track.close();
		interpolate_from(std::valarray<double>(1.0, 1),
				std::list< std::valarray<double> >(1,
					list_to_valarray(quantities[AGE])),
				std::list< std::valarray<double> >(1,
					list_to_valarray(quantities[RSTAR])),
				std::list< std::valarray<double> >(1,
					list_to_valarray(quantities[ICONV])),
				std::list< std::valarray<double> >(1,
					list_to_valarray(quantities[IRAD])),
				std::list< std::valarray<double> >(1,
					list_to_valarray(quantities[MRAD])),
				std::list< std::valarray<double> >(1,
					list_to_valarray(quantities[RRAD])),

				smoothing[RSTAR], smoothing[ICONV], smoothing[IRAD],
				smoothing[MRAD], smoothing[RRAD],

				nodes[RSTAR], nodes[ICONV], nodes[IRAD], nodes[MRAD],
				nodes[RRAD],
				
				(column_numbers[LSTAR]==-1 ?
				 std::list< std::valarray<double> >() :
				 std::list< std::valarray<double> >(1,
					 list_to_valarray(quantities[LSTAR]))),
				smoothing[LSTAR], nodes[LSTAR],
				std::numeric_limits<double>::max(), 0);
	}
};
