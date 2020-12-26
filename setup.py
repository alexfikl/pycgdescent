import sys
from setuptools import setup

from pybind11.setup_helpers import Pybind11Extension, build_ext

# TODO: add support for conditionally using blas; currently does not build on
# readthedocs with it on

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
                )
            ],
        cmd_class={"build_ext": build_ext},
        )
