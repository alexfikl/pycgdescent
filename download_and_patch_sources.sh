#!/usr/bin/env bash
set -Eeuo pipefail

pkgname='CG_DESCENT-C'
pkgver='6.8'
archive="${pkgname}-${pkgver}.tar.gz"
url="https://users.clas.ufl.edu/hager/papers/CG/Archive/${archive}"

basedir=$(pwd)
builddir="${basedir}/_build"
patchdir="${basedir}/patches"

# {{ get original sources

mkdir -p "${builddir}"
pushd "${builddir}"

echo -e "\033[1;32mbuilddir: $(pwd)\033[0m"

if [ ! -f "${builddir}/${archive}" ]; then
    echo -e "\033[1;32mDownloading...\033[0m"
    wget -q -N "${url}"
fi

echo -e "\033[1;32mExtracting '${archive}'...\033[0m"
tar xvf "${archive}"

# }}}

# {{{ apply patches

declare -a patches=(
    '0001-add-extern-c.patch'
    '0002-add-header-guards.patch'
    '0003-add-func-typedefs.patch'
    '0004-add-user-pointer-to-functions.patch'
)

pushd "${pkgname}-${pkgver}"
pwd
for patch in "${patches[@]}"
do
    echo -e "\033[1;32mApplying '${patchdir}/${patch}\033[0m'"
    patch -p1 -i "${patchdir}/${patch}"
done

# }}}

# {{{ copy sources

echo -e "\033[1;32mCopying patched sources...\033[0m"
for filename in cg_user.h cg_blas.h cg_descent.h cg_descent.c
do
    cp "${filename}" "${basedir}/src"
done

popd
popd
# rm -rf ${builddir}

# }}}
