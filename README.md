# Simplified Betts-Miller convection scheme

[![Build and test](https://github.com/brian-rose/climlab-sbm-convection/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/brian-rose/climlab-sbm-convection/actions/workflows/build-and-test.yml)


Brian Rose, University at Albany

_This is a work in progress!_

## About

This is a stand-alone Python wrapper for the Simplified Betts-Miller convection scheme described by [Frierson (2007), J. Atmos. Sci. 64, doi:10.1175/JAS3935.1](https://doi.org/10.1175/JAS3935.1).

Thanks to Dargan Frierson for sharing the original Fortran code.


## Building from source

### Build environment

Here are instructions to create a build environment (including Fortran compiler) with conda/mamba

Starting from the root of the `climlab-sbm-convection` repo *(example for Apple M1 machine, see `./ci/` for other environment files)*:
```
mamba env create --file ./ci/requirements-macos-arm64.yml
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

For now, see the notebook in this repo.

An appropriate runtime environment for this notebook can be found [here](https://github.com/brian-rose/ClimateLaboratoryBook/blob/main/environment.yml).
