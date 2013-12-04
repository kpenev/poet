#ifndef __IO_UTIL_H
#define __IO_UTIL_H

#include <iostream>
#include <valarray>
#include <list>

///Exctracts columns from a white space separated stream.
template<typename COLNUM_STRUCTURE>
std::valarray< std::list<double> > parse_columns(
		///The stream which should contain only columns, and the columns
		///with number listed in column_numbers must be interpretable as
		///doubles.
		std::istream &is,

		///The column numbers to extract. Should provide operator[], size and
		///iterators.
		const COLNUM_STRUCTURE &column_numbers);

template<class COLNUM_STRUCTURE>
std::valarray< std::list<double> > parse_columns(std::istream &is,
		const COLNUM_STRUCTURE &column_numbers)
{
	int last_column=*max_element(column_numbers.begin(),
			column_numbers.end());
	std::string line;
	std::valarray< std::list<double> > result(column_numbers.size());
	unsigned line_number=0;
	for(std::getline(is, line); !is.eof(); std::getline(is, line)) {
		++line_number;
		if(!is.good()) {
			std::ostringstream msg;
			msg << "Failed to read line " << line_number
				<< " of the data section of an input stream in "
				<< "parse_columns.";
			throw Error::IO(msg.str());
		}
		if(line[0]=='#') continue;
		std::istringstream line_stream(line);
		for(int column_number=0; column_number<=last_column;
				++column_number) {
			std::ostringstream msg;
			msg << "Failed to parse column " << column_number
				<< " on line " << line_number << " of the data section of an"
				" input stream in parse_columns.";
			if(!line_stream.good()) throw Error::IO(msg.str());
			typename COLNUM_STRUCTURE::const_iterator
				colnum_iter=column_numbers.begin();
			for(size_t result_index=0; result_index<column_numbers.size();
					++result_index, ++colnum_iter)
				if(*colnum_iter==column_number) {
					double value;
					line_stream >> value;
					result[result_index].push_back(value);
					break;
				}
			if(colnum_iter==column_numbers.end()) {
				std::string word;
				line_stream >> word;
			}
			if(line_stream.bad()) throw Error::IO(msg.str());
		}
	}
	return result;
}

#endif
