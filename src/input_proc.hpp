#pragma once
#include "core_mesh.hpp"
#include "pin_mesh.hpp"

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