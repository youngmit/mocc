# MOCC root SConstruct
pugixml_include = Dir('#lib/pugixml/src')

# env = Environment(CXX = 'clang++', CXXFLAGS="-std=c++11 -m64",
#                   LINKFLAGS = "-m64",
#                   CPPPATH=[pugixml_include])
env = Environment(CXXFLAGS="-std=c++11",
                  CPPPATH=[pugixml_include])
Export('env')

VariantDir('build', 'src', duplicate=0)
VariantDir('build/lib', 'lib')
lib_objs = SConscript(['build/lib/SConscript'])

Export('lib_objs')

SConscript(['build/SConscript'])