# MOCC root SConstruct
import os

AddOption(
    '--debug-build',
    dest='debug-build',
    action='store_true',
    help='build with debug symbols',
    default=False)

AddOption(
    '--profile-build',
    dest='profile',
    action='store_true',
    help='build with profiling information',
    default=False)


pugixml_include = Dir('#lib/pugixml/src')

cxx = 'clang++'

cxxflags = "-std=c++11 -Wall"
linkflags = "-lboost_regex"

if GetOption('debug-build'):
    cxxflags += " -g"
if GetOption('profile'):
    cxxflags += " -pg"
    linkflags += " -pg"
    cxx = "g++"

env = Environment(CXX=cxx,
                  CXXFLAGS=cxxflags,
                  LINKFLAGS=linkflags,
                  CPPPATH=[pugixml_include])


env['ENV']['TERM'] = os.environ['TERM']

Export('env')

VariantDir('build', 'src', duplicate=0)
VariantDir('build/lib', 'lib')
lib_objs = SConscript(['build/lib/SConscript'])

Export('lib_objs')

SConscript(['build/SConscript'])
