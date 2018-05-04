#ifndef COMMAND_H
#define COMMAND_H


#include <iostream>
#include <deque>
#include <string>
#include "stdio.h" 

class Command
{
public:
	Command(const char *name)
	{
		commandline.assign(name);
		prefetched.assign(0,name);
		const size_t start( commandline.find_first_of("([{\t\n") );
		const size_t end( commandline.find_first_of(")]}\n", start+1) );
		if(start < commandline.length())
		{
			head.assign( commandline.substr(0, start) );
			if( end > (start + 1) )
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
		
		
		/*size_t found = commandline.find_first_of("([{\t\n");
		size_t end = commandline.find_last_of(")]}\n");
		head = commandline.substr(0, found);
		args = commandline.substr(found + 1, end - found - 1);*/
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

	void next(FILE* input);
	void push_back(std::string const& newcommand);
	bool assertCommand(FILE* input, const char *name, int val);
	bool assertCommand(FILE* input, const char *name, long val);
	bool assertCommand(FILE* input, const char *name, double val);

	
	bool getCommand(FILE* input, const char *name, int &val);
	bool getCommand(FILE* input, const char *name, long &val);
	bool getCommand(FILE* input, const char *name, double &val);



protected:


private:
	std::string commandline;
	std::deque<std::string> prefetched;
	std::string head;
	std::string args;
};


#endif
