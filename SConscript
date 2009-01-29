# -*- python -*-

tools = Split("""
qt4
profilerdds
doxygen
qtt_qtconfig
boost_program_options
""")

env = Environment(tools = ['default'] + tools)

qt4modules = ['QtCore','QtGui']
env.EnableQt4Modules(qt4modules)

sources = Split("""
main.cpp
ScopeTSReader.cpp
ProfilerScope.cpp
""")

headers = Split("""
ScopeTSReader.h
ProfilerScope.h
""")

html = env.Apidocs(sources + headers, DOXYFILE_FILE = "Doxyfile")

profilerscope = env.Program('profilerscope', sources)

Default(profilerscope, html)

