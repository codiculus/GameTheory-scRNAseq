from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize


setup(
  name = "combinations",
  ext_modules = cythonize("combinations.pyx")
  )