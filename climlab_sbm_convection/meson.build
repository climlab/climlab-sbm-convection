python_sources = [
  '__init__.py',
]

py.install_sources(
  python_sources,
  pure: false,
  subdir: 'climlab_sbm_convection'
)

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

mysource = custom_target('_simplified_betts_miller.c',
  input : ['climlab_betts_miller.f90'],  
  output : ['_simplified_betts_millermodule.c'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', 
             '-m', '_simplified_betts_miller', '--lower', 
             '--build-dir', '@OUTDIR@']
)

inc_np = include_directories(incdir_numpy, incdir_f2py)

py.extension_module('_simplified_betts_miller',
  ['climlab_betts_miller.f90', mysource],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  subdir: 'climlab_sbm_convection'
)

subdir('tests')