# -*- python -*-

tools = Split("""
qt4
qtt_scopeplot
qtt_knob
doxygen
fftw
""")

env = Environment(tools = ['default'] + tools)

qt4modules = ['QtCore','QtGui']
env.EnableQt4Modules(qt4modules)

# This will create ui_AScope.h
env.Uic4(['AScope.ui',])

sources = Split("""
AScope.cpp
PlotInfo.cpp
""")

headers = Split("""
AScope.h
PlotInfo.h
""")

html = env.Apidocs(sources + headers, DOXYFILE_DICT={'PROJECT_NAME':'ascope', 'PROJECT_NUMBER':'1.0'})

ascope = env.Library('ascope', sources)

Default(ascope, html)

tooldir = env.Dir('.').srcnode().abspath    # this directory

def ascope(env):
    env.AppendUnique(CPPPATH = [tooldir])
    env.AppendLibrary('ascope')
    env.AppendDoxref('ascope')
    env.Require(tools)

Export('ascope')



