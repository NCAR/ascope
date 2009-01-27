# -*- python -*-

tools = Split("""
qt4
profilerdds
doxygen
""")

env = Environment(tools = ['default'] + tools)

qt4modules = ['QtCore',]
env.EnableQt4Modules(qt4modules)

sources = Split("""
main.cpp
""")

headers = Split("""
""")

html = env.Apidocs(sources + headers, DOXYFILE_FILE = "Doxyfile")

profilerscope = env.Program('profilerscope', sources)

Default(profilerscope, html)

