#pragma once
#include "coremesh.hpp"
#include "pinmesh.hpp"

#include <map>
#include <memory>

class InputProc {
public:
	InputProc(const char* filename);
	void GenerateMesh();
private:
	std::map<int, SP_PinMesh> m_pinMeshes;
};