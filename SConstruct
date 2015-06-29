# MOCC root SConstruct
pugixml_include = Dir('#lib/pugixml/src')
boost_include = Dir('/usr/local/include/')

env = Environment(CXXFLAGS="-std=c++11",
                  LINKFLAGS="-L/usr/local/lib /usr/local/lib/libboost_regex.a",
                  CPPPATH=[pugixml_include, boost_include])
#env = Environment(CXXFLAGS="-std=c++11",
#                  CPPPATH=[pugixml_include])
Export('env')

VariantDir('build', 'src', duplicate=0)
VariantDir('build/lib', 'lib')
lib_objs = SConscript(['build/lib/SConscript'])

Export('lib_objs')

SConscript(['build/SConscript'])
