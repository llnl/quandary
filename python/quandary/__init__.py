"""Quandary Python interface for quantum optimal control."""

import warnings

from .quandary import *
from . import new

warnings.warn(
    "The legacy Quandary interface is deprecated and will be removed in a "
    "future version. Please migrate to the new interface: "
    "'from quandary.new import Quandary, QuandaryConfig'",
    DeprecationWarning,
    stacklevel=2
)
