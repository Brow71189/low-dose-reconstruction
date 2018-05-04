#include "stdio.h"
#include "string.h"
#include <sstream>
#include <string>

#include "command.h"
#include "globals.h"
#include "fast_cache.h"

void Command::next(FILE* input)
{
	//printf("reading input ... ");
	if(prefetched.empty())
	{
		do
		{
			commandline.clear();
			///These one prevents a strange bug when running without cctable
			///it is most likely some memory corruption in update_subsums_from_black()
			///it only shows up after the second iteration.
			commandline.reserve(128);
			do
			{
				int ch(getc(input));
				if(ch != EOF)
				{	commandline.push_back(ch); }
				else
				{
					printf("reached end of input stream\n");
					commandline.clear();
					commandline = "Run(0)\n";
					break;	
				}
			} while (commandline.back() != '\n');
		//skip empty lines and comments	
		} while ((commandline.length() <= 1) || (commandline[0]=='#') );
	}
	else
	{
		commandline = prefetched.front();
		prefetched.pop_front();
	}
	const size_t start( commandline.find_first_of("([{\t\n") );
	const size_t end( commandline.find_first_of(")]}\n", start+1) );
	if( (start != std::string::npos) && (start < commandline.length()) ) //either one should suffice
	{
		head.assign( commandline.substr(0, start) );
		if( end != std::string::npos && end > (start + 1) )
		{
			args.assign( commandline.substr(start + 1, end - start - 1) );
		}
		else//either an empty argument, unpared brackets or whatever
		{
			args.clear();
		}
	}
	else //should be impossible
	{
		head.assign(commandline);
		args.clear();
	}
	return;
}

void Command::push_back(std::string const& newcommand)
{
	prefetched.push_back(newcommand);
}

bool Command::getCommand(FILE* input, const char *name, int &val)
{
	next(input);
	if( (strcmp(name, head.c_str()) != 0) ||
		 (sscanf(args.c_str(), "%d", &val) != 1) )
	{
		printf( "%s(%s) != %s(int)\n", head.c_str(), args.c_str(), name);
		return false;
	}
	return true;     
}

bool Command::getCommand(FILE* input, const char *name, long &val)
{
	next(input);
	if( (strcmp(name, head.c_str()) != 0) ||
		 (sscanf(args.c_str(), "%ld", &val) != 1) )
	{
		printf( "%s(%s) != %s(long)\n", head.c_str(), args.c_str(), name);
		return false;
	}
	return true;     
}
	
bool Command::getCommand(FILE* input, const char *name, double &val)
{
	next(input);
	if( (strcmp(name, head.c_str()) != 0) ||
		 (sscanf(args.c_str(), "%lf", &val) != 1) )
	{
		printf( "%s(%s) != %s(float)\n", head.c_str(), args.c_str(), name);
		return false;
	}
	return true;     
}

bool Command::assertCommand(FILE* input, const char *name, int val)
{
	int readout = -1;
	next(input);
	if( (strcmp(name, head.c_str()) != 0) ||
		(sscanf(args.c_str(), "%d", &readout) != 1) ||
		 (readout != val) )
	{
		printf( "%s(%s) != %s[%d] ?\n", head.c_str(), args.c_str(), name, val);
		return false;
	}
	return true;
}

bool Command::assertCommand(FILE* input, const char *name, long val)
{
	long readout = -1;
	next(input);
	if( (strcmp(name, head.c_str()) != 0) ||
		(sscanf(args.c_str(), "%ld", &readout) != 1) ||
		 (readout != val) )
	{
		printf( "%s(%s) != %s[%ld] ?\n", head.c_str(), args.c_str(), name, val);
		return false;
	}
	return true;
}


bool Command::assertCommand(FILE* input, const char *name, double val)
{
	double readout = -1;
	next(input);
	if( (strcmp(name, head.c_str()) != 0) ||
		(sscanf(args.c_str(), "%lf", &readout) != 1) ||
		 ((int)readout != (int)val) )
	{
		printf( "%s(%s) != %s[%lf] ?\n", head.c_str(), args.c_str(), name, val);
		return false;
	}
	return true;
}

