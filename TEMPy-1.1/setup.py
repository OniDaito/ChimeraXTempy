import sys, os,shutil
from distutils.core import setup

if os.path.exists('build')==True:
	print "build exists"
	shutil.rmtree('./build')

try:
    import numpy
except ImportError:
    sys.stderr.write('numpy is not installed, you can find it at: '
                     'http://numpy.scipy.org\n')
    sys.exit()

#control version
if [int(dgt) for dgt in numpy.__version__.split('.')[:2]] < [1, 6]:
    sys.stderr.write('TEMPy requires numpy v1.6 or later, you can find it at: '
                     'http://numpy.scipy.org\n')
    sys.exit()

#biopython>=1.51/1.4?
try:
    import Bio
except ImportError:
    sys.stderr.write('Biopython is not installed, you can find it at: '
                     'http://biopython.org/wiki/Main_Page\n')
    sys.exit()

if [int(dgt) for dgt in Bio.__version__.split('.')[:2]] < [1, 5]:
    sys.stderr.write('TEMPy requires Biopython v1.5 or later, you can find it at: '
                     'http://biopython.org/wiki/Main_Page\n')
    sys.exit()

#Scipy>=???
try:
    import scipy
except ImportError:
    sys.stderr.write('Scipy is not installed, you can find it at: '
                     'http://www.scipy.org/\n')
    sys.exit()

if [int(dgt) for dgt in scipy.__version__.split('.')[:2]] < [0, 1]:
    sys.stderr.write('TEMPy requires Scipy v0.1 or later, you can find it at: '
                     'http://www.scipy.org/\n')
    sys.exit()

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 5):
    print "TEMPy requires Python 2.5 or better. Python %d.%d detected" % \
        sys.version_info[:2]
    print "Please upgrade your version of Python."
    sys.exit(-1)


setup(
    name='TEMPy',
    version='1.0',
    author='Maya Topf, Daven Vasishtan, Arun Prasad Pandurangan, Irene Farabella, Agnel-Praveen Joseph, Harpal Sahota',
    author_email='tempy-help.cryst.bbk.ac.uk',
    packages=['TEMPy'],
    url='http://tempy.ismb.lon.ac.uk/',
#    license='LICENSE.txt',
    description='TEMPy: a Python Library for Assessment of 3D Electron Microscopy Density Fits',
#    long_description=open('README.txt').read(),
    package_dir = {'TEMPy':'src/TEMPy'},
    package_data = {'TEMPy': ['tempy_data/*']},
    requires=['NumPy (>=1.6)',
        "Scipy (>= 0.1)",
        "Biopython (>= 1.5)",
    ],
)
