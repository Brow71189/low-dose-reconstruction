#ifndef COMMAND_HPP
#define COMMAND_HPP


#include <iostream>
#include <deque>
#include <string>
#include "stdio.h" 

//#define BUFFER_LEN 40

class Command
{
public:
	Command(const char *name)
	{
		commandline.reserve(128);
		head.reserve(128);
		args.reserve(16);
		
		commandline.assign(name);
		const size_t start( commandline.find_first_of("([{\t\n") );
		const size_t end( commandline.find_first_of(")]}\n", start+1) );
		if( (start != std::string::npos) && (start < commandline.length()) )
		{
			head.assign( commandline.substr(0, start) );
			if( (end != std::string::npos) && (end > (start + 1)) )
			{
				args.assign( commandline.substr(start + 1, end - start - 1) );
			}
			else//either an empty argument, unpared brackets or whatever
			{
				args.clear();
			}
		}
		else
		{
			head.assign(commandline);
			args.clear();
		}
		
	};
	~Command(void){};
	const std::string& getHead(void)// The bare head
	{
		return head;
	};
	const std::string& getArgs(void)// The naked argumnets
	{
		return args;
	};

	void next();
	void push_back(std::string const& newcommand);
	bool assertCommand(const char *name, int val);
	bool assertCommand(const char *name, double val);

	bool getCommand(const char *name, int &val);
	bool getCommand(const char *name, double &val);



protected:


private:
	std::string commandline;
	std::deque<std::string> prefetched;
	std::string head;
	std::string args;
};


#endif
