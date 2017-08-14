# Tempy for ChimeraX

Welcome to the ChimeraXTempy bundle plugin. This project combines tempy with ChimeraX for scoring EM fitting.

## History & Version
Currently at version 0.0.5

## Building

This plugin for ChimeraX is pure-python though it does install BioPython as part of the setup.

### Building and testing for development

#### ChimeraX approach

This is the accepted approach starting from 0.0.2. With ChimeraX setup, change the paths in the Makefile to the location of your chimera install. With that in place, you can execute the following:

    make install

It should run fine, downloading biopython and installing automatically. You can then run:

    ./ChimeraX --start Tempy --cmd 'open ~/Projects/ChimeraXTempy/test/final1_mdcg.pdb; open ~/Projects/ChimeraXTempy/test/1akeA_10A.mrc;'

#### Virtualenv approach

Deprecated in light of the above! Setup a virtualenv environment with the lib/python3.6 directory set to ChimeraX's lib/python3.6 directory. This way, the build_install command can use pip to install the bundle into the right place. Secondly, alter the two paths at the top of the build_install.sh script.

Inside the virtualenv setup with ChimeraX

    ./build_install.sh

To load the UI and some test models, use something like

    ./ChimeraX --start Tempy --cmd 'open ~/Projects/ChimeraXTempy/test/final1_mdcg.pdb; open ~/Projects/ChimeraXTempy/test/1akeA_10A.mrc;'

### Bundle to distribute

This is a little easier. Simple execute the following:

    make wheel

This creates a python wheel bundle in the dist directory. You can then use this with ChimeraX to install the plugin.

## How to use this plugin

### Windowed tool

Go to the following menu

    General -> Tempy

The window should appear in your general set of tools, usually to the right of the screen.

To score with SCCC you need to select one map, one model and make sure you have set the *rigid file* dialog. Assuming that model 1 is a model and model 2 is a volume map:

    sel #1 #2 

The *rigid file* button will launch a file dialog so you can select the correct rigid file for scoring. You can then press the *sccc* button. If you want to change the resolution of the map (in Angstroms) you can do this from the *Res.* textbox. *sigma* can also be set from the neighbouring text box.

The SMOC score can accept multiple models. To use this feature, you can perform multiple selects:

    sel #1 #2 #3

Each model is scored against a single map, and plotted in the tempy dialog box. Depending on the size of your screen the graph might not appear. You may need to stretch the tempy window if it is docked.

The SMOC score accepts the *Res.*, *sigma* and *window* parameters. Once you have selected at least one model and one map, press the *SMOC* button.

Scores and messages appear in the *log* screen.

The NMI score can be calculated by selecting two models, a model and a map or two maps. Clicking on the the NMI score button will print the resulting score in the *log* screen. The variables for the nmi score include *Res.* *Res.2* (the resolution of the first and second selected items respectively), *c1* and *c2* (the contour values for the first and second selected items respectively).

### Command-line tools
Inside ChimeraX with the command line, you can compute the sccc with the following command

    sccc (#1) (#2) ~/Projects/ChimeraXTempy/test/rigid_RF.txt

Where #1 is the pdb model, #2 is the map. The SMOC is similar:

    smoc (#1,2) (#3)

You can pass in multiple models as the first parameter in order to compare scores.

The NMI score works in the same way

    nmi (#1) (#2) 10 10 6 6

The first two parameters are the selected items. The third and fourth are the resolutions, and the last two, the contours.

## Testing

The test directory included contains a python test and some sample model and map files for verifying the scoring functions. These can be run using the command:

    python test_cmd.py

This assumes you have the ChimeraX python modules in your path, as well as TEMPY. The easiest way to to do this is to use VirtualEnv or similar.

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

### BioPython

All the fine folks who make BioPython.

[https://github.com/biopython/](https://github.com/biopython/)
