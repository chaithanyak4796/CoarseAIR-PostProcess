#ifndef LOGGER_H
#define LOGGER_H

// The functions here just print all the inputs in one line.


//static std :: stringstream Log;
static void Write()
{
         std ::cout<<std::endl;
	//printf("%s \n",Log.str());
	//Log.str(std::string());
}


template <typename T, typename... Types> 
static void Write(T var1, Types... var)
{
        std :: cout<<var1<<" ";
	//Log<<var1<" ";
	Write(var...);
	
}

#endif /* LOGGER_H */
