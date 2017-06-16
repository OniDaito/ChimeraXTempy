#!/bin/bash
# Assuming inside a virtualenv that reflects the ChimeraX lib python site-packages dir
HOME=/home/oni
CHIMERAX=/opt/UCSF/chimerax
rm -rf ./build
rm -rf ./dist
rm -rf $HOME/.cache/ChimeraX/0.1/toolshed && rm -rf $CHIMERAX/lib/python3.6/site-packages/ChimeraX_Tempy-0.0.1.dist-info && rm -rf $CHIMERAX/lib/python3.6/site-packages/chimerax/tempy
make
pip install dist/ChimeraX_Tempy-0.0.1-py3-none-any.whl
