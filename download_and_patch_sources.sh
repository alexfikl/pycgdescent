#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2020-2024 Alexandru Fikl <alexfikl@gmail.com>
#
# SPDX-License-Identifier: MIT

set -Eeuo pipefail

# {{{ description

pkgname='SuiteOPT'
pkgver='3.0.2'
archive="${pkgname}-${pkgver}.tar_.gz"
url="https://people.clas.ufl.edu/hager/files/${archive}"

basedir=$(pwd)
builddir="${basedir}/build"
patchdir="${basedir}/patches"

# }}}

mkdir -p "${builddir}"
pushd "${builddir}"

# {{{ get original sources

echo -e "\033[1;32mbuilddir: $(pwd)\033[0m"

if [ ! -f "${builddir}/${archive}" ]; then
  echo -e "\033[1;32mDownloading...\033[0m"
  curl -L -O "${url}"
fi

echo -e "\033[1;32mExtracting '${archive}'...\033[0m"
tar xvf "${archive}"

# }}}

# {{{ patch sources

pushd "${pkgname}"

declare -a patches=(
  '0001-feat-guard-BLAS-defines.patch'
  '0002-feat-extern-C-in-cgdescent-header.patch'
  '0003-feat-add-user-callback.patch'
)

for patch in "${patches[@]}"; do
  echo -e "\033[1;32mApplying '${patchdir}/${patch}\033[0m'"
  patch -p1 -i "${patchdir}/${patch}"
done

popd

# }}}

# {{{ copy sources

echo -e '\033[1;32mCopying patched sources...\033[0m'
mkdir -p ${basedir}/src/wrapper

declare -a files=(
  'SuiteOPT/CGDESCENT/Include/cg_descent.h'
  'SuiteOPT/CGDESCENT/Source/cg_default.c'
  'SuiteOPT/CGDESCENT/Source/cg_descent.c'
  'SuiteOPT/CGDESCENT/Source/cg_print.c'
  'SuiteOPT/CGDESCENT/Source/cg_util.c'
  'SuiteOPT/SSM/Include/SSM.h'
  'SuiteOPT/SSM/Source/SSM.c'
  'SuiteOPT/SSM/Source/SSMdiagopt.c'
  'SuiteOPT/SSM/Source/SSMmult.c'
  'SuiteOPT/SSM/Source/SSMprint.c'
  'SuiteOPT/SSM/Source/SSMrefine.c'
  'SuiteOPT/SSM/Source/SSMtridiag.c'
  'SuiteOPT/SuiteOPTconfig/sopt.h'
  'SuiteOPT/SuiteOPTconfig/sopt.c'
)

for filename in "${files[@]}"; do
  cp "${filename}" "${basedir}/src/wrapper"
done

# }}}

popd
# rm -rf ${builddir}

# vim:set ts=2 sts=2 sw=2 et:
