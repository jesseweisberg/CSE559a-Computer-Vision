
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: users.cpp                                                                            //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"

Users::Users()
:
std::list<User>()
{
	//empty
}

void Users::save(std::string filename) const
{
	BinaryFileWriter file(filename);
	save(file);
}

void Users::save(BinaryFileWriter& file) const
{
	file.write(getSize());
	for(int i=0; i<getSize(); i++) {
		(*this)[i].save(file);
	}
	std::cout << "Saved users to '" << file.getFilename() << "'" << std::endl;

}

void Users::load(std::string filename)
{
	BinaryFileReader file(filename);
	load(file);

}

int Users::getSize() const
{
	return int(std::list<User>::size());
}

void Users::load(BinaryFileReader& file)
{	
	int size=file.readInt();
	for(int i=0; i<size; i++) {
		User new_user;
		new_user.load(file);
		addUser(new_user);
	}
	//std::cout << "Loaded users from '" << file.getFilename() << "'" << std::endl;

}

void Users::addUser(const User& user)
{
	push_back(user);
	map_string[user.getName()]=&back();
	map_int[int(size()-1)]=&back();
}

const User& Users::operator[](int index) const
{

	std::map<int, const User*>::const_iterator i;
	if ((i=map_int.find(index))==map_int.end()) {
		throw Error("Index '%%' does not correspond to a known user", index); 
	}
	return (*(*i).second);
}

const User& Users::operator[](std::string name) const
{
	std::map<std::string, const User*>::const_iterator i;
	if ((i=map_string.find(name))==map_string.end()) {
		throw Error("Name '%%' does not correspond to a known user", name); 
	}
	return (*(*i).second);

}

void Users::sort()
{
	std::list<User>::sort();
	int n=0;
	for(std::list<User>::const_iterator it=begin(); it!=end(); it++, n++) {
		map_string[(*it).getName()]=&(*it);
		const User* user=&(*it);
		map_int[n]=user;

	}

}



