project('climlab_sbm_convection', 'c')

py = import('python').find_installation(pure: false)

py.extension_module(
    'our_first_module',
    'our_first_module.c',
    install: true,
)

py3.extension_module('_simplified_betts_miller',
  flapack_module,
  c_args: [Wno_empty_body],
  link_args: version_link_args,
  dependencies: [lapack_dep, blas_dep, fortranobject_dep],
  install: true,
  subdir: 'scipy/linalg'
)