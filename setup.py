import sys
from setuptools import setup

from pybind11.setup_helpers import Pybind11Extension, build_ext

sources = [
        "src/cg_descent.c",
        "src/cg_descent_wrap.cpp"
        ]

setup(
        ext_modules=[
            Pybind11Extension("_cg_descent",
                sources=sources,
                # NOTE: to enable use of std::optional
                cxx_std=17,
                # FIXME: more robust way to find blas
                extra_link_args=["-lm", "-lblas", "-lpthread"]),
            ],
        cmd_class={"build_ext": build_ext},
        )
