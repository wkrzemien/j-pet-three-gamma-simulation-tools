#ifndef Message_h
#define Message_h
#include <iostream>
#include <string>

class Message
{
	public:
		static bool Check(const bool v, std::string msg);
		static bool Error(const bool v, std::string msg);
};

#endif
