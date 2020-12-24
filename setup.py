import sys
from setuptools import setup

from pybind11.setup_helpers import Pybind11Extension, build_ext

sources = [
        "src/cg_descent.c",
        "src/cg_descent_wrap.cpp"
        ]

setup(
        ext_modules=[
            Pybind11Extension("_cg_descent", sources,
                # NOTE: to enable use of std::optional
                extra_compile_args=["-std=c++17"]),
            ],
        cmd_class={"build_ext": build_ext},
        )
