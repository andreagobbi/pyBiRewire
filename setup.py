from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy as np

NAME = "BiRewire"
VERSION = "0.2"
DESCR = "A Cython wrapper of BiRewire"
URL = "http://bioconductor.org/packages/devel/bioc/html/BiRewire.html"
REQUIRES = ['numpy', 'cython']

AUTHOR = "Andrea Gobbi"
EMAIL = "gobbi.andrea@mail.com"

LICENSE = "GPL V3"

SRC_DIR = "BiRewire"
PACKAGES = [SRC_DIR]

ext_1 = Extension(SRC_DIR + ".wrapped",
                  [SRC_DIR + "/lib/BiRewire.c", SRC_DIR + "/wrapped.pyx"],
                  libraries=[],
                  include_dirs=[np.get_include()])


EXTENSIONS = [ext_1]

if __name__ == "__main__":
    setup(install_requires=REQUIRES,
          packages=PACKAGES,
          zip_safe=False,
          name=NAME,
          version=VERSION,
          description=DESCR,
          author=AUTHOR,
          author_email=EMAIL,
          url=URL,
          license=LICENSE,
          cmdclass={"build_ext": build_ext},
          ext_modules=EXTENSIONS
          )
