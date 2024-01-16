# SPDX-FileCopyrightText: 2020-2022 Alexandru Fikl <alexfikl@gmail.com>
#
# SPDX-License-Identifier: MIT

import os

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

# {{{ enable blas

CG_DISABLE_BLAS = bool(int(os.environ.get("CG_DISABLE_BLAS", 0)))
if CG_DISABLE_BLAS:
    extra_link_args = []
    defines = [("CG_DISABLE_BLAS", None)]
else:
    # FIXME: better way to find / choose blas library?
    extra_link_args = ["-lblas", "-lpthread"]
    defines = []

# }}}

sources = ["src/wrapper/cg_descent.c", "src/wrapper/cg_descent_wrap.cpp"]

setup(
    ext_package="pycgdescent",
    ext_modules=[
        Pybind11Extension(
            "_cg_descent",
            sources=sources,
            # NOTE: to enable use of std::optional
            cxx_std=17,
            extra_link_args=extra_link_args,
            define_macros=defines,
        )
    ],
    cmdclass={"build_ext": build_ext},
)
