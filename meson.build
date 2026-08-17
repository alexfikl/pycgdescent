project(
    'pycgdescent',
    'c',
    'cpp',
    license: 'GPL-2.0-or-later',
    meson_version: '>=1.3.0',
    default_options: ['buildtype=release', 'warning_level=2', 'cpp_std=c++17'],
    subproject_dir: 'src/subprojects',
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

nanobind_dep = dependency('nanobind')
if get_option('use-blas')
    blas_dep = dependency(['openblas', 'OpenBLAS', 'cblas'], required: false)
    if not blas_dep.found()
        warning('No BLAS implementation found.')
        add_project_arguments('-DCG_DISABLE_BLAS', language: 'c')
    endif
else
    blas_dep = []
    add_project_arguments('-DCG_DISABLE_BLAS', language: 'c')
endif

py.extension_module(
    '_cg_descent',
    ['src/wrapper/cg_descent.c', 'src/wrapper/cg_descent_wrap.cpp'],
    dependencies: [nanobind_dep, blas_dep],
    c_args: ['-Wno-maybe-uninitialized'],
    limited_api: '3.12',
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
