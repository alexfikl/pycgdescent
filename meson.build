project('pycgdescent', 'c', 'cpp',
  license: 'GPL-2.0-or-later',
  meson_version: '>=1.1.0',
  default_options : [
    'buildtype=debug',
    'warning_level=3',
    'c_std=c18',
    'cpp_std=c++17',
    'use-blas=true',
  ],
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

pybind11_dep = dependency('pybind11', version: '>=2.11.0')
if get_option('use-blas')
  blas_dep = dependency(['openblas', 'OpenBLAS', 'cblas'], required: false)
  if not blas_dep.found()
    add_project_arguments('-DSUITEOPT_DISABLE_BLAS', language: 'c')
  endif
else
  blas_dep = []
  add_project_arguments('-DSUITEOPT_DISABLE_BLAS', language: 'c')
endif

add_project_arguments([
  '-D_DEFAULT_SOURCE',
  '-DSUITEOPT_BLAS_UNDERSCORE',
  # '-DDLONG',
  ], language: 'c')

cg_descent_sources = [
  # sources
  'src/wrapper/SSM.c',
  'src/wrapper/SSMdiagopt.c',
  'src/wrapper/SSMmult.c',
  'src/wrapper/SSMprint.c',
  'src/wrapper/SSMrefine.c',
  'src/wrapper/SSMtridiag.c',
  'src/wrapper/cg_default.c',
  'src/wrapper/cg_descent.c',
  'src/wrapper/cg_print.c',
  'src/wrapper/cg_util.c',
  'src/wrapper/sopt.c',
  # wrapper
  'src/wrapper/cg_descent_wrap.cpp'
]

py.extension_module(
  '_cg_descent', cg_descent_sources,
  dependencies: [pybind11_dep, blas_dep],
  c_args: ['-DDLONG'],
  install: true,
  subdir: 'pycgdescent')

py.install_sources([
  'src/pycgdescent/__init__.py',
  'src/pycgdescent/py.typed',
  'src/pycgdescent/_cg_descent.pyi',
  ],
  subdir: 'pycgdescent')
