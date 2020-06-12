"""
Setup for epi_inference package
"""

import os
from setuptools import setup, find_packages

requires=[ 'pyomo', 'pyyaml', 'jsondiff' ]

setup(name="epi_inference",
    version='1.0.dev0',
    maintainer='Carl D. Laird',
    maintainer_email='cdlaird@sandia.gov',
    platforms = ["any"],
    python_requires='>=3.6',
    description = 'Tools to perform inference for epidemiological models',
    long_description = """
This package was developed to estimate disease transmission parameters
associated with epidemiology models.  The package includes several
inference approaches.

This package relies on the Pyomo modeling software.  Pyomo and other
optimization solvers need to be installed to use these inference methods.
""",
    long_description_content_type='text/markdown',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Mathematics'],
      packages=find_packages(),
      keywords=['utility'],
      install_requires=requires,
      entry_points="""
        [console_scripts]
        epiinf = epi_inference.epiinf:main
      """
      )

