#pragma once
#include <iostream>
#include <stdlib.h>
#include <string>

using std::cout;
using std::endl;
using std::string;

static void Error(const char* msg){
	cout << "ERROR: " << msg << endl;
	exit(EXIT_FAILURE);
}
static void Warn(const char* msg){
	cout << "WARNING: " << msg << endl;
}