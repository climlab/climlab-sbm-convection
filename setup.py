from os.path import join, abspath

VERSION = '0.1'

def configuration(parent_package='', top_path=None):
    global config
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.fcompiler import get_default_fcompiler, CompilerNotFound

    build = True
    try:
        # figure out which compiler we're going to use
        compiler = get_default_fcompiler()
        # set some fortran compiler-dependent flags
        f90flags = []
        if compiler == 'gnu95':
            f90flags.append('-fno-range-check')
            f90flags.append('-ffree-form')
        elif compiler == 'intel' or compiler == 'intelem':
            f90flags.append('-132')
        #  Set aggressive optimization level
        f90flags.append('-O3')
        #  Suppress all compiler warnings (avoid huge CI log files)
        f90flags.append('-w')
    except CompilerNotFound:
        print('No Fortran compiler found, not building the SimplifiedBettsMiller module!')
        build = False

    config = Configuration(package_name='_simplified_betts_miller', parent_name=parent_package, top_path=top_path)
    if build:
        config.add_extension(name='_simplified_betts_miller',
                             sources=[gen_source],
                             extra_f90_compile_args=f90flags,
                             f2py_options=['--quiet'],
                             )
    return config

def gen_source(ext, build_dir):
    thispath = config.local_path
    sourcelist = []
    sourcelist.append(join(thispath,'climlab_betts_miller.f90'))
    #sourcelist.append(join(thispath,'convect.f'))
    #sourcelist.append(join(thispath,'Driver.f90'))
    try:
        config.have_f90c()
        return sourcelist
    except:
        print('No Fortran 90 compiler found, not building SimplifiedBettsMiller extension!')
        return None

def setup_package():
    __version__ = VERSION
    metadata = dict(
          name = '_simplified_betts_miller',
          version = __version__,
          setup_requires=['numpy'],
          license='MIT',
          )
    run_build = True
    from setuptools import setup
    if run_build:
        from numpy.distutils.core import setup
        metadata['configuration'] = configuration
    setup(**metadata)

if __name__ == '__main__':
    setup_package()
