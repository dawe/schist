from . import inference
from . import tools as tl
from . import plotting as pl
from . import io


__author__ = ', '.join([
    'Davide Cittaro',
    'Leonardo Morelli',
])
__email__ = ', '.join([
    'cittaro.davide@hsr.it',
    'l.morelli2@studenti.unisr.it',
])

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
