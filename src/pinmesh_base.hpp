#pragma once

#include "globalconfig.hpp"

class PinMesh{
public:
	int id(){
		return m_id;
	}
protected:
	int m_id;
	mocc::float_t m_pitchX;
	mocc::float_t m_pitchY;
};