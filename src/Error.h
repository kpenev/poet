/**\addtogroup Utilities_group
 * @{
 * 
 * \file
 *
 * \brief Defines the exception hierarchy used by this code.
 *
 */

#ifndef __ERROR_H
#define __ERROR_H

#include <iostream>
#include <exception>
#include <string>

///Isolates all exceptions.
namespace Error {
	///\brief The base class of all exceptions.
	///
	///Supports what() which should return the type of error and a
	///separate get_message() which gives details about the error.
	class General : public std::exception {
	private:
		std::string message;
	public:
		///Create an exception.
		General(const std::string &error_message="") : 
			message(error_message) {}
		
		///The type of error.
		virtual const char *what() const throw() {return "general error";}

		///Details about what caused the error.
		const std::string &get_message() {return message;}

		///Cleanup if necessary.
		virtual ~General() throw() {}
	};

	///%Error detected by the ALGLIB library.
	class ALGLIB : public General {
	public:
		///Create an exception for an error detected inside ALGLIB.
		ALGLIB(const std::string &error_message="") : 
			General(error_message) {}
		
		///Reports "ALGLIB error" as the error type.
		virtual const char *what() const throw()
		{return "ALGLIB error";}
	};

	///Any runtime error.
	class Runtime : public General {
	public:
		///Create a runtime exception.
		Runtime(const std::string &error_message="") :
			General(error_message) {}

		///Reports "Run time error" as the eror type.
		virtual const char *what() const throw()
		{return "Run time error";}
	};

	///Function arguments do not satisfy some requirement.
	class BadFunctionArguments : public Runtime {
	public:
		///Create bad function arguments exception.
		BadFunctionArguments(const std::string &error_message="") :
			Runtime(error_message) {
		}

		///Reports "Bad function arguments" as the erorr type.
		virtual const char *what() const throw()
		{return "Bad function arguments";}
	};

	///Exception indicating unrecognized or unsuitable stellar zone.
	class BadStellarZone : public BadFunctionArguments {
	public:
		///Create a bad stellar zone exception.
		BadStellarZone(const std::string &error_message="") :
			BadFunctionArguments(error_message) {}

		///Reports "Invalid stellar zone" as the error type.
		virtual const char *what() const throw()
		{return "Invalid stellar zone";}
	};

	///Exception indicating that a file or a directory was not found.
	class PathNotFound : public Runtime {
	private:
		///Whether the problem was with a directory and not a file.
		bool directory;
	public:
		///Create a missing file or directory exception.
		PathNotFound(const std::string &message, 
			     const std::string &filename="", 
			     bool isdir=false) :
			Runtime(filename+", "+message), directory(isdir) {}

		///Reports "File/Directory not found" as the error type.
		virtual const char *what() const throw() {return (directory ?
				"Directory not found." : 
				"File not found.");}
	};

	///Input/Output exception.
	class IO : public Runtime {
	public:
		///Create an Input/Ouput exception.
		IO(const std::string &filename="", 
				const std::string &error_message="") :
			Runtime("Error reading "+filename+", "
					+error_message) {}
		///Reports "I/O error" as the error type
		virtual const char *what() const throw() {return "I/O error";}
	};

	///%Error related to parsing the command line.
	class CommandLine : public Runtime {
	public:
		///Create command line parsing exception.
		CommandLine(const std::string &error_message="") :
			Runtime(error_message) {}

		///Returns "Error parsing the command line" as the error type
		virtual const char *what() const throw()
		{return "Error parsing the command line";}
	};
};
#endif

///@}
