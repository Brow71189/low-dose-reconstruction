
#include "globals.hpp"
#include "command.hpp"
#include <assert.h>

void Command::next()
{
	//printf("reading input ... ");
	try
	{
		if(prefetched.empty())
		{
			do
			{
				commandline.clear();
				do
				{
					int ch( getc(stdin) );
					if( (ch != EOF) ) // && (ch != '\0') 
					{	commandline.push_back(ch); }
					else
					{
						printf("reached end of input stream\n");
						fflush(stdout);
						fprintf(logfile,"reached end of input stream / connection lost\n defaulting to \"Run(0)\\n\"\n");
						fflush(logfile);
						commandline.assign("Run(0)\n");
						//break;	
					}
				} while (commandline.back() != '\n');
			//skip empty lines and comments		
			} while ( (commandline.length() <= 1) || (commandline[0]=='#') );
		}
		else
		{
			//assert(false); //we dont expect any prefetching
			commandline = prefetched.front();
			prefetched.pop_front();
		}
		
	}
	catch(...)
	{
		commandline.assign("ReadError(stdin)\n");
	}
	fprintf(logfile,"%s\t",commandline.substr(0,commandline.length()-1).c_str());fflush(logfile);
	const size_t start( commandline.find_first_of("([{\t\n") );
	const size_t end( commandline.find_first_of(")]}\n", start+1) );
	if( (start != std::string::npos) && (start < commandline.length()) ) //either one should suffice
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
	else //should be impossible
	{
		printf("ERROR: a command without defined ending was encountered");
		fflush(stdout);
		head.assign(commandline);
		args.clear();
	}
	fprintf(logfile,"%s[%s]\n",head.c_str(),args.c_str());fflush(logfile);
	return;
}

void Command::push_back(std::string const& newcommand)
{
	prefetched.push_back(newcommand);
}

bool Command::getCommand(const char *name, int &val)
{
	next();
	if( (strcmp(name, head.c_str()) != 0) ||
		 (sscanf(args.c_str(), "%d", &val) != 1) )
	{
		printf( "%s(%s) != %s[int_]\n", head.c_str(), args.c_str(), name);
		fflush(stdout);
		return false;
	}
	return true;     
}
	
bool Command::getCommand(const char *name, double &val)
{
	next();
	if( (strcmp(name, head.c_str()) != 0) ||
		 (sscanf(args.c_str(), "%lf", &val) != 1) )
	{
		printf( "%s(%s) != %s[float_]\n", head.c_str(), args.c_str(), name);
		fflush(stdout);
		return false;
	}
	return true;     
}

bool Command::assertCommand(const char *name, int val)
{
	int readout = -1;
	next();
	if( (strcmp(name, head.c_str()) != 0) ||
		(sscanf(args.c_str(), "%d", &readout) != 1) ||
		 (readout != val) )
	{
		printf( "%s(%s) != %s[%d] ?\n", head.c_str(), args.c_str(), name, val);
		fflush(stdout);
		return false;
	}
	return true;
}

bool Command::assertCommand(const char *name, double val)
{
	double readout = -1;
	next();
	if( (strcmp(name, head.c_str()) != 0) ||
		(sscanf(args.c_str(), "%lf", &readout) != 1) ||
		 ((int)readout != (int)val) )
	{
		printf( "%s(%s) != %s[%lf] ?\n", head.c_str(), args.c_str(), name, val);
		fflush(stdout);
		return false;
	}
	return true;
}

