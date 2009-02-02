# -*- python -*-

tools = Split("""
qt4
qtt_qtconfig
qtt_scopeplot
qtt_knob
profilerdds
doxygen
fftw
boost_program_options
""")

env = Environment(tools = ['default'] + tools)

qt4modules = ['QtCore','QtGui']
env.EnableQt4Modules(qt4modules)

# This will create ui_ProfilerScope.h
env.Uic4(['ProfilerScope.ui',])

sources = Split("""
main.cpp
ScopeTSReader.cpp
ProfilerScope.cpp
PlotInfo.cpp
""")

headers = Split("""
ScopeTSReader.h
ProfilerScope.h
PlotInfo.h
""")

html = env.Apidocs(sources + headers, DOXYFILE_FILE = "Doxyfile")

profilerscope = env.Program('profilerscope', sources)

Default(profilerscope, html)



