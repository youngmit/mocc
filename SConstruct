# MOCC root SConstruct
import os

AddOption(
    '--debug-build',
    dest='debug-build',
    action='store_true',
    help='build with debug symbols',
    default=False)

pugixml_include = Dir('#lib/pugixml/src')

cxx = 'clang++'

if GetOption('debug-build'):
    env = Environment(CXX=cxx,
                  CXXFLAGS="-g -std=c++11",
                  LINKFLAGS="-lboost_regex",
                  CPPPATH=[pugixml_include])
else:
    env = Environment(CXXFLAGS="-std=c++11",
                  LINKFLAGS="-lboost_regex",
                  CPPPATH=[pugixml_include])

env['ENV']['TERM'] = os.environ['TERM']

Export('env')

VariantDir('build', 'src', duplicate=0)
VariantDir('build/lib', 'lib')
lib_objs = SConscript(['build/lib/SConscript'])

Export('lib_objs')

SConscript(['build/SConscript'])
