import os
import sys

from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

# {{{ enable blas

CG_DISABLE_BLAS = bool(int(os.environ.get("CG_DISABLE_BLAS", 0)))
if CG_DISABLE_BLAS:
    extra_link_args = []
    defines = [("CG_DISABLE_BLAS", None)]
else:
    # FIXME: better way to find / choose blas library?
    extra_link_args = ["-blas", "-lpthread"]
    defines = []

# }}}

sources = [
        "src/cg_descent.c",
        "src/cg_descent_wrap.cpp"
        ]

setup(
        ext_package="pycgdescent",
        ext_modules=[
            Pybind11Extension("_cg_descent",
                sources=sources,
                # NOTE: to enable use of std::optional
                cxx_std=17,
                extra_link_args=extra_link_args,
                define_macros=defines,
                )
            ],
        cmdclass={"build_ext": build_ext},
        )
