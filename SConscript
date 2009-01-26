# -*- python -*-

tools = Split("""
qt4
profilerdds
""")

env = Environment(tools = ['default'] + tools)

qt4modules = ['QtCore',]
env.EnableQt4Modules(qt4modules)

sources = Split("""
main.cpp
""")

p = env.Program('profilerscope', sources)
Default(p)

