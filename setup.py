from os import cpu_count
from setuptools import setup
from pathlib import Path
from setuptools import setup, find_packages
import versioneer

_path = Path('requirements.txt')
with _path.open() as requirements:
    requires = [l.strip() for l in requirements]


    
try:
    import graph_tool.all as gt
except ImportError:
    raise ImportError(
            """Please install the graph-tool library either visiting

            https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions

            or by conda: `conda install -c conda-forge graph-tool`
            """
    )    

setup(name='schist',
      description='Nested Stochastic Block Model applied to single cell data.',
      long_description='Partitioning of neighborhood graphs in Scanpy using Nested Stochastic Block Model.',
      url='http://github.com/dawe/scNSBM',
      author='Davide Cittaro',
      author_email='cittaro.davide@gmail.com',
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
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',          
      ],
      python_requires='>=3.4',
      zip_safe=False,
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      )
