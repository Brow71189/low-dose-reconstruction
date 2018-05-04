#ifndef SPARSE_ARRAY_HEADER
#define SPARSE_ARRAY_HEADER

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "globals.h"


//pretty much the same bad performance for either approach
//#define USE_VECTOR_AND_SET


#ifdef USE_VECTOR_AND_SET
class Sparse_Array
{
public:
	Sparse_Array(unsigned int new_default = 0)
	{
		usual = new_default;
	};
	~Sparse_Array(void){};
	void put_val(unsigned int position, unsigned int value)
	{
		if(value == usual)
		{	return;}
		value -= ( (value>usual) ? 1 : 0 );
		if(data.size() < (value+1) )
		{	data.resize(value+1);}
		data[value].insert(position);
	}
	unsigned int get_val(unsigned int position)
	{
		for(unsigned int i(0); i<data.size(); ++i)
		{
			if(data[i].count(position))
			return (i>usual) ? (i+1) : i ;
		}
		return usual;
	}
	unsigned int count(void)
	{
		unsigned int total(0);
		for(unsigned int i(0); i < data.size();++i)
		{
			total += data[i].size();
		}
		return total;
	}

protected:

private:
	unsigned int usual;
	std::vector<std::unordered_set<unsigned int> > data;
};
#else
class Sparse_Array
{
public:
	Sparse_Array(unsigned int new_default = 0)
	{
		usual = new_default;
	};
	~Sparse_Array(void){};
	void put_val(unsigned int position, unsigned int value)
	{
		if(value == usual)
		{	return;}
		data[position] = value;
	}
	unsigned int get_val(unsigned int position)
	{
		std::unordered_map<unsigned int, unsigned int>::iterator
		it( data.find(position) );
		return (it != data.end() ) ? it->second : usual;
	}
	unsigned int count(void)
	{
		return data.size();
	}

protected:

private:
	unsigned int usual;
	std::unordered_map<unsigned int, unsigned int> data;
};

#endif

#endif
