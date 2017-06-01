# Tempy for ChimeraX

Welcome to the ChimeraXTempy bundle plugin. This project combines tempy with ChimeraX for scoring EM fitting.

## History & Version
Currently at version 0.0.1

## Building

### Building and testing for development

Setup a virtualenv environment with the lib/python3.6 directory set to ChimeraX's lib/python3.6 directory. This way, the build_install command can use pip to install the bundle into the right place. Secondly, alter the two paths at the top of the build_install.sh script.

Inside the virtualenv setup with ChimeraX

    ./build_install.sh

To run, something like

    ./ChimeraX --start Tempy --cmd 'open ~/Projects/tempy/data/tempy/1akeA_10A.mrc; open ~/Projects/tempy/data/tempy/1akeA.pdb'

### Bundle to distribute

This is a little easier. Simple execute the following:

    make

This creates a python wheel bundle in the dist directory. You can then use this with ChimeraX to install the plugin.


## Tempy

Taken from [http://tempy.ismb.lon.ac.uk http://tempy.ismb.lon.ac.uk] 

"TEMPy is an object-oriented Python library designed to help the user in the analysis of structures of macromolecular assemblies, especially in the context of 3D electron microscopy density maps. It is designed with a set of functionalities that assess the goodness-of-fit between a given atomic model and a density map or between two maps using a variety of different scoring functions.It can also generate various ensembles of alternative fits, which has been shown to be useful in assessing a model fit and find other models with good density fit. In the future, TEMPy will also include a suite of functions for density ﬁtting."


## Credits

### ChimeraX

From [https://www.cgl.ucsf.edu/chimerax/ https://www.cgl.ucsf.edu/chimerax/]

ChimeraX is the work of many individuals:

- Thomas Ferrin - principal investigator
- Conrad Huang - project leader
- Eric Pettersen - atomic structure analysis, sequence-structure linkage
- Tom Goddard - density maps, hierarchical models, rendering
- Greg Couch - infrastructure, I/O
- Elaine Meng - user documentation, web content
- Scooter Morris - project management 

... and numerous collaborators, to whom we are grateful for thoughtful suggestions and time spent evaluating features. 

### tempy
Copyright © 2015 Birkbeck College University of London.
Authors: Maya Topf, Daven Vasishtan, Arun Prasad Pandurangan, Irene Farabella, Agnel-Praveen Joseph, Harpal Sahota

    Farabella, I., Vasishtan, D., Joseph, A.P., Pandurangan, A.P., Sahota, H. & Topf, M.
    TEMPy: a Python library for assessment of three-dimensional electron microscopy 
    density fits. (2015). J. Appl. Cryst. 48, 1314-1323, link.
