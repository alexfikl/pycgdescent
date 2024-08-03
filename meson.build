project(
    'pycgdescent',
    'c',
    'cpp',
    license: 'GPL-2.0-or-later',
    meson_version: '>=1.1.0',
    default_options: ['buildtype=release', 'warning_level=2', 'cpp_std=c++17'],
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

pybind11_dep = dependency('pybind11', version: '>=2.11.0')
if get_option('use-blas')
    blas_dep = dependency(['openblas', 'OpenBLAS', 'cblas'], required: false)
    if not blas_dep.found()
        add_project_arguments('-DCG_DISABLE_BLAS', language: 'c')
    endif
else
    blas_dep = []
    add_project_arguments('-DCG_DISABLE_BLAS', language: 'c')
endif

py.extension_module(
    '_cg_descent',
    ['src/wrapper/cg_descent.c', 'src/wrapper/cg_descent_wrap.cpp'],
    dependencies: [pybind11_dep, blas_dep],
    install: true,
    subdir: 'pycgdescent',
)

py.install_sources(
    [
        'src/pycgdescent/__init__.py',
        'src/pycgdescent/py.typed',
        'src/pycgdescent/_cg_descent.pyi',
    ],
    subdir: 'pycgdescent',
)
