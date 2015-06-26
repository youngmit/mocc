#pragma once
#include "coremesh.hpp"
#include "pinmesh.hpp"

#include <map>
#include <memory>

namespace mocc{

class InputProc {
public:
	InputProc(const char* filename);
	void GenerateMesh();
private:
	std::map<int, SP_PinMesh> m_pinMeshes;
};

};