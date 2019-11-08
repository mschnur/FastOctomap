#!/bin/bash

THIS_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

pushd "$THIS_SCRIPT_DIR"

cat rgbd_dataset_freiburg1_desk2.tar.gz.parta* >rgbd_dataset_freiburg1_desk2.tar.gz
tar xvzf rgbd_dataset_freiburg1_desk2.tar.gz

popd