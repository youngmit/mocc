# MOCC root SConstruct
env = Environment(tools = ['mingw'])

Export('env')

VariantDir('build', 'src', duplicate=0)
VariantDir('build/lib', 'lib')
lib_objs = SConscript(['build/lib/SConscript'])

Export('lib_objs')

SConscript(['build/SConscript'])