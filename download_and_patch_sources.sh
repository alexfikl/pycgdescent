#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2020-2022 Alexandru Fikl <alexfikl@gmail.com>
#
# SPDX-License-Identifier: MIT

set -Eeuo pipefail

# {{{ description

pkgname='CG_DESCENT-C'
pkgver='6.8'
archive="${pkgname}-${pkgver}.tar_.gz"
url="https://people.clas.ufl.edu/hager/files/${archive}"

basedir=$(pwd)
builddir="${basedir}/build"
patchdir="${basedir}/patches"

# }}}

# {{{ get original sources

mkdir -p "${builddir}"
pushd "${builddir}"

echo -e "\033[1;32mbuilddir: $(pwd)\033[0m"

if [ ! -f "${builddir}/${archive}" ]; then
  echo -e "\033[1;32mDownloading...\033[0m"
  curl -L -O "${url}"
fi

echo -e "\033[1;32mExtracting '${archive}'...\033[0m"
tar xvf "${archive}"

# }}}

# {{{ apply patches

declare -a patches=(
  '0000-add-blas-compile-flag.patch'
  '0001-add-extern-c.patch'
  '0002-add-header-guards.patch'
  '0003-add-func-typedefs.patch'
  '0004-add-user-pointer-to-functions.patch'
  '0005-add-iteration-callback.patch'
  '0006-add-step-size-limit.patch'
  '0007-cg_evaluate-initialize-df.patch'
)

pushd "${pkgname}-${pkgver}"
pwd
for patch in "${patches[@]}"; do
  echo -e "\033[1;32mApplying '${patchdir}/${patch}\033[0m'"
  patch -p1 -i "${patchdir}/${patch}"
done

# }}}

# {{{ copy sources

echo -e "\033[1;32mCopying patched sources...\033[0m"
mkdir -p ${basedir}/src/wrapper

for filename in cg_user.h cg_blas.h cg_descent.h cg_descent.c; do
  cp "${filename}" "${basedir}/src/wrapper"
done

popd
popd
# rm -rf ${builddir}

# }}}

# vim:set ts=2 sts=2 sw=2 et:
