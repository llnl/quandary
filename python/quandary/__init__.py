"""Quandary Python interface for quantum optimal control."""

import warnings

from .quandary import *
from . import new

warnings.warn(
    "The legacy Quandary class interface is deprecated and will be removed in a "
    "future version. If you are using 'from quandary.new import *', you can safely "
    "ignore this warning.",
    DeprecationWarning,
    stacklevel=2
)
