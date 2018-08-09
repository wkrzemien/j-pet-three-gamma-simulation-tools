#include "Message.h"

bool Message::Check(const bool v, std::string msg)
{
	if(!v)
		std::cout<<msg<<std::endl;
	return v;
}

bool Message::Error(const bool v, std::string msg)
{
	return !Check(v, "[ERROR] " + msg);
}
