#!/bin/bash

ver=$(cat find_interstit.sh | grep version | head -1 | gawk -F "=" '{print $2}')

echo "ver= "$ver

files=$(find -follow -type f | grep -v make_package.sh | grep -v find_interstitials.tgz | grep -v archive)

tar  chvzf package_find_interstitials_v$ver.tgz $files

