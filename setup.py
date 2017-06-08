# vim: set expandtab ts=4 sw=4:

# import distutils.core
# distutils.core.DEBUG = True
# TODO: remove distlib monkey patch when the wheel package
# implements PEP 426's pydist.json
from distlib import metadata
metadata.METADATA_FILENAME = "metadata.json"
from setuptools import setup, Extension
import os, os.path, sys
	
description = """Tempy."""
	
# Assume Python executable is in ROOT/bin/python
# and make include directory be ROOT/include
root = os.path.dirname(os.path.dirname(sys.executable))
inc_dir = os.path.join(root, "include")
lib_dir = os.path.join(root, "lib")

if sys.platform == "darwin":
    # Tested with macOS 10.12
    libraries = []
    compiler_flags = ["-std=c++11", "-stdlib=libc++"]
    env = "Environment :: MacOS X :: Aqua",
    op_sys = "Operating System :: MacOS :: MacOS X"
elif sys.platform == "win32":
    # Tested with Cygwin
    libraries = ["libatomstruct"]
    compiler_flags = []
    env = "Environment :: Win32 (MS Windows)"
    op_sys = "Operating System :: Microsoft :: Windows :: Windows 10"
else:
    # Presumably Linux
    # Tested with Ubuntu 16.04 LTS running in
    #   a singularity container on CentOS 7.3
    libraries = []
    compiler_flags = ["-std=c++11"]
    env = "Environment :: X11 Applications"
    op_sys = "Operating System :: POSIX :: Linux"
	
# If platform-specific, list C++ source files
#ext_sources = ["src/_sample.cpp"]
# If pure Python, set as empty list
ext_sources = []
	
if ext_sources:
    # platform-specific
    ext_mods = [Extension("chimerax._tempy",
        define_macros=[("MAJOR_VERSION", 0),
                       ("MINOR_VERSION", 1)],
                        extra_compile_args=compiler_flags,
                        include_dirs=[inc_dir],
                        library_dirs=[lib_dir],
                        libraries=libraries,
                        sources=ext_sources)]
    environments = [env, op_sys]
else:
    # pure Python
    ext_mods = []
    environments = [
        "Environment :: MacOS X :: Aqua",
        "Environment :: Win32 (MS Windows)",
        "Environment :: X11 Applications",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows :: Windows 10",
        "Operating System :: POSIX :: Linux",
    ]

setup(
    name="ChimeraX_Tempy",
    version="0.0.1",  # PEP 440, should match Development Status below
    description="The Tempy plugin for ChimeraX",  # one line synopsis
    long_description=description,  # see above
    author="Birkbeck",
    author_email="oni@section9.co.uk",
    url="http://tempy.ismb.lon.ac.uk/",
    python_requires=">= 3.5",
    package_dir={
        "chimerax.tempy": "src",
        "TEMPy" : "TEMPy",
        # Wonder if we should put TEMPy in here?
    },
    packages=[
        "chimerax.tempy",
        "TEMPy",
    ],
    ext_modules=ext_mods,
    install_requires=[
        # list dependences on non-standard Python packages incl. ChimeraX ones
        # Question - should we package TEMPy with this or load externally here?
        # For now I'm including it as TEMPy isn't in pip
        "ChimeraX-Core",
        "numpy",
        "scipy",
        "biopython",
    ],
    classifiers=[
        # From https://pypi.python.org/pypi?%3Aaction=list_classifiers
        # and our own ChimeraX classifiers.
        "Development Status :: 2 - Pre-Alpha",  # TODO: update as appropriate
        "Framework :: ChimeraX",
        "Intended Audience :: Science/Research",
        "License :: Free for non-commercial use",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "ChimeraX :: Bundle :: General :: 1,1 :: chimerax.tempy :: chimerax.test :: ",
        "ChimeraX :: Tool :: Tempy :: General :: Show Tempy",
        "ChimeraX :: Command :: sccc :: General :: Calculate the SCCC with Tempy",
    ] + environments,
)
