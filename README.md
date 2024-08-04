# Simplified Betts-Miller convection scheme

[![Build and test](https://github.com/climlab/climlab-sbm-convection/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/climlab/climlab-sbm-convection/actions/workflows/build-and-test.yml)


Brian Rose, University at Albany

## About

This is a stand-alone Python wrapper for the Simplified Betts-Miller convection scheme described by [Frierson (2007), J. Atmos. Sci. 64, doi:10.1175/JAS3935.1](https://doi.org/10.1175/JAS3935.1).

Thanks to Dargan Frierson for sharing the original Fortran code.

The primary use-case is to serve as a moist convection driver option
for [climlab](https://climlab.readthedocs.io/), but it can also be used 
as a stand-alone model. This is a lightweight wrapper that emulates the 
Fortran interface as closely as possible. 
Dargan Frierson's original Fortran code is bundled here in the `src` directory for reference.

## Installing

Pre-built binaries for many platforms are available from [conda-forge](https://conda-forge.org).

To install in the current environment:
```
conda install climlab-sbm-convection --channel conda-forge
```
or create a self-contained environment:
```
conda create --name my_env python=3.11 climlab-sbm-convection --channel conda-forge
conda activate my_env
```

See below for instructions on how to build from source.

## Building from source

### Build environment

Here are instructions to create a build environment (including Fortran compiler) with conda/mamba

Starting from the root of the `climlab-sbm-convection` repo *(example for Apple M1 machine, see `./ci/` for other environment files)*:
```
mamba env create --file ./ci/requirements-macos-arm64.yml
conda activate sbm_build_env
```

Or, to specify the Python version, you can do
```
mamba create --name sbm_build_env python=3.11 --channel conda-forge
mamba env update --file ./ci/requirements-macos-arm64.yml
conda activate sbm_build_env
```

### Building and installing into the Python environment

From the root of the repository, do this:
```
python -m pip install . --no-deps -vv
```

### Running tests

To run tests, do this from any directory other than the climlab-emanuel-convection repo:
```
pytest -v --pyargs climlab_sbm_convection
```

##  Example usage

For now, see the notebooks in the `docs` directory of this repo.

An appropriate runtime environment for these notebooks is defined in the file `docs/example_env.yml`. To create this environment and launch the notebook server, navigate to the `docs` directory and do this:

```
conda env create -f example_env.yml
conda activate sbm_example_env
jupyter lab
```

## Version history

- Version 0.2 (released 8/2/2024) implements support for a non-uniform relative humidity parameter. The subroutine `betts_miller` now expects array input for relative humidity. The call signature for this subroutine has changed.
- Version 0.1 is the first public release. 