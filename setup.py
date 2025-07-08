from os import cpu_count
from setuptools import setup
from pathlib import Path
from setuptools import setup, find_packages

_path = Path('requirements.txt')
with _path.open() as requirements:
    requires = [l.strip() for l in requirements]

setup(name='schist',
      description='Nested Stochastic Block Model applied to single cell data.',
      long_description='Partitioning of neighborhood graphs in Scanpy using Nested Stochastic Block Model.',
      url='http://github.com/dawe/schist',
      license='BSD 3',
      packages=find_packages(),
      install_requires=requires,
      dependency_links=[
          'https://git.skewed.de/count0/graph-tool/tree/master'
      ],
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 3.8',          
          'Programming Language :: Python :: 3.10',          
          'Programming Language :: Python :: 3.12',          
      ],
      )
