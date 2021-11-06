#include<iostream>
#include <cstdlib>
#include<fstream>
#include<math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include "Main.h"
#include "Global.h"

using namespace std;

int Find_length(const std :: string& file_name, int rskip)
{
	// This function simply finds the length of each file.
	
	int len = 0;
	std :: string line;
	ifstream flen;
	flen.open(file_name.c_str());
	
	if(!flen)
	{
		cout<<" Error opening file : "<<file_name<<endl;
		exit(0);
	}
	
	while(std::getline(flen,line))
	{
		len ++;
	}
	
	flen.close();
	
	len = len - rskip;  // Number of rows to skip
	
	return len;
}
