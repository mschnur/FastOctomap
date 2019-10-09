#!/bin/bash

THIS_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

function build_octomap_lib() { 
    local src_dir="$1"
	local build_dir="$2"
	local install_dir="$3"
	local common_cmake_options="-DBUILD_OCTOVIS_SUBPROJECT=OFF -DBUILD_DYNAMICETD3D_SUBPROJECT=OFF -DOCTOVIS_QT5=OFF"
	
	
	mkdir -p "${build_dir}"
	
	pushd "${build_dir}"
	cmake -DCMAKE_INSTALL_PREFIX="${install_dir}" ${common_cmake_options} "${src_dir}"
    make -j$(nproc) install
	make -j$(nproc) test
	popd
};

ROOT_BUILD_DIR="${THIS_SCRIPT_DIR}/build"
ORIGINAL_OCTOMAP_BUILD_DIR="${ROOT_BUILD_DIR}/original_octomap"
MOD_OCTOMAP_BUILD_DIR="${ROOT_BUILD_DIR}/octomap"

ROOT_INSTALL_DIR="${THIS_SCRIPT_DIR}/deploy"
ORIGINAL_OCTOMAP_INSTALL_DIR="${ROOT_INSTALL_DIR}/original_octomap"
MOD_OCTOMAP_INSTALL_DIR="${ROOT_INSTALL_DIR}/octomap"

build_octomap_lib "$THIS_SCRIPT_DIR/original_octomap" "$ORIGINAL_OCTOMAP_BUILD_DIR" "$ORIGINAL_OCTOMAP_INSTALL_DIR"
build_octomap_lib "$THIS_SCRIPT_DIR/octomap" "$MOD_OCTOMAP_BUILD_DIR" "$MOD_OCTOMAP_INSTALL_DIR"
