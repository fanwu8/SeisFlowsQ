#!/bin/sh

rm -rf output*
rm -rf scratch

SFROOT="$(dirname ""$(dirname "$PWD")"")"
export PYTHONPATH=$SFROOT
python $SFROOT"/scripts/sfrun"