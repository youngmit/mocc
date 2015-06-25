#pragma once
// This provides a number of utility classes for doing timing. For now there is
// just a simple class providing tic/toc functionality.

class Timer{
public:
	Timer();
	tic();
	float toc();	
};