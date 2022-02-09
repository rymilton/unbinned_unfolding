import sys

try:
    from skbuild import setup
except ImportError:
    print('Please update pip, you need pip 10 or greater,\n'
          ' or you need to install the PEP 518 requirements in pyproject.toml yourself', file=sys.stderr)
    raise

setup(
    name="RooUnfold",
    cmake_args=['-DBUILD_DOC=off'],
    version='2.0.0',
    url='https://gitlab.cern.ch/RooUnfold/RooUnfold',
    author='Giordon Stark',
    author_email='gstark@cern.ch',
    description='pip-installable package of RooUnfold',
    packages=['RooUnfold'],
    package_dir={'': 'src'}
)
