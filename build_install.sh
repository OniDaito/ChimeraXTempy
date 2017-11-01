#!/bin/bash
# Assuming inside a virtualenv that reflects the ChimeraX lib python site-packages dir
VERSION=0.3.0
HOME=/home/oni
CHIMERAX=/opt/UCSF/chimerax
rm -rf ./build
rm -rf ./dist
rm -rf $HOME/.cache/ChimeraX/0.1/toolshed && rm -rf $CHIMERAX/lib/python3.6/site-packages/ChimeraX_Tempy-$VERSION.dist-info && rm -rf $CHIMERAX/lib/python3.6/site-packages/chimerax/tempy
make
pip install dist/ChimeraX_Tempy-$VERSION-py3-none-any.whl
