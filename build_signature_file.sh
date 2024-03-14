#  Run this script to re-generate the signature file _simplified_betts_miller.pyf
f2py --overwrite-signature climlab_betts_miller.f90 -m _simplified_betts_miller -h _simplified_betts_miller.pyf
