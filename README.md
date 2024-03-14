# Simplified Betts-Miller convection scheme

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
mamba create --name sbm_build_env python=3.10 --channel conda-forge
mamba env update --file ./ci/requirements-macos-arm64.yml
conda activate sbm_build_env
```

### Building with f2py

From the root of the repository, do this:
```
f2py -c -m _simplified_betts_miller climlab_betts_miller.f90
```

This will create the shared object `_simplified_betts_miller.cpython-*.so` that can be imported in a Python session.

_This is a work in progress, and eventually there will be better packaging so this manual build step will not be needed._

##  Example usage

For now, see the notebook in this repo.

An appropriate runtime environment for this notebook can be found [here](https://github.com/brian-rose/ClimateLaboratoryBook/blob/main/environment.yml).
