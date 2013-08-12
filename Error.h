#ifndef __ERROR_H
#define __ERROR_H
#include <iostream>
#include <exception>
#include <string>

namespace Error {
	class General : public std::exception {
	private:
		std::string message;
	public:
		General(const std::string &error_message="") : 
			message(error_message) {}
		virtual const char *what() const throw() {return "general error";}
		const std::string &get_message() {return message;}
		virtual ~General() throw() {}
	};

	class ALGLIB : public General {
	public:
		ALGLIB(const std::string &error_message="") : 
			General(error_message) {}
		virtual const char *what() const throw() {return "ALGLIB error";}
	};

	class Runtime : public General {
	public:
		Runtime(const std::string &error_message="") :
			General(error_message) {}
		virtual const char *what() const throw() {return "Run time error";}
	};

	class BadFunctionArguments : public Runtime {
	public:
		BadFunctionArguments(const std::string &error_message="") :
			Runtime(error_message) {
		}
		virtual const char *what() const throw()
		{return "Bad function arguments";}
	};

	class BadStellarZone : public BadFunctionArguments {
	public:
		BadStellarZone(const std::string &error_message="") :
			BadFunctionArguments(error_message) {}
		virtual const char *what() const throw()
		{return "Invalid stellar zone";}
	};

	class PathNotFound : public Runtime {
	private:
		bool directory;
	public:
		PathNotFound(const std::string &message, 
			     const std::string &filename="", 
			     bool isdir=false) :
			Runtime(filename+", "+message) {}
		virtual const char *what() const throw() {return (directory ? 
				"Directory not found." : 
				"File not found.");}
	};

	class IO : public Runtime {
	public:
		IO(const std::string &filename="", 
				const std::string &error_message="") :
			Runtime("Error reading "+filename+", "
					+error_message) {}
		virtual const char *what() const throw() {return "I/O error";}
	};

	class CommandLine : public Runtime {
	public:
		CommandLine(const std::string &error_message="") :
			Runtime(error_message) {}
		virtual const char *what() const throw()
		{return "Error parsing the command line";}
	};
};
#endif
