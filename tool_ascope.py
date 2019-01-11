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

env['DOXYFILE_DICT'].update({'PROJECT_NAME':'Ascope'})
html = env.Apidocs(sources + headers)

ascope = env.Library('ascope', sources)

Default(ascope)

tooldir = env.Dir('.').srcnode().abspath    # this directory

def ascope(env):
    env.AppendUnique(CPPPATH = [tooldir])
    env.AppendLibrary('ascope')
    env.AppendDoxref('ascope')
    env.Require(tools)

Export('ascope')



