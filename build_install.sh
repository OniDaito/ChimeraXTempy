#!/bin/bash
 rm -rf /home/oni/.cache/ChimeraX/0.1/toolshed && rm -rf ~/Projects/chimerax/lib/python3.6/site-packages/ChimeraX_Tempy-0.0.1.dist-info && rm -rf ~/Projects/chimerax/lib/python3.6/site-packages/chimerax/tempy
 make
 pip install dist/ChimeraX_Tempy-0.0.1-py3-none-any.whl
