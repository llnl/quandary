"""Python subclasses of C++ nanobind-bound structs.

Each class in this module wraps a C++ configuration struct and adds:

- **Numpy conversion**: numpy arrays, integers, and floats are automatically
  converted to Python builtins before being passed to the C++ layer, so
  ``setup.nlevels = np.array([2, 3])`` works seamlessly.
- **Improved error messages**: ``TypeError`` messages include the field name,
  the rejected value, and a hint when a negative integer is assigned to an
  unsigned field.

See the C++ docstrings on each field (accessible via ``help()``) for
type information and allowed values.
"""

import numpy as np

from .._quandary_impl import Setup as _CppSetup
from .._quandary_impl import (
    InitialConditionSettings as _CppInitialConditionSettings,
    OptimTargetSettings as _CppOptimTargetSettings,
    ControlParameterizationSettings as _CppControlParameterizationSettings,
    ControlInitializationSettings as _CppControlInitializationSettings,
)


def _to_builtin(value):
    """Recursively convert numpy types to Python builtins.

    Nanobind's type-casters expect native Python types (``list``, ``int``,
    ``float``), so this function converts numpy arrays and scalars before
    they reach the C++ property setters.
    """
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, list):
        return [_to_builtin(v) for v in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def _fmt_val(v):
    """Return a repr of *v* truncated to 80 characters for error messages."""
    r = repr(v).replace('\n', ' ')
    return r[:77] + '...' if len(r) > 80 else r


def _make_repr(cpp_cls):
    """Build a ``__repr__`` that introspects read-write fields on *cpp_cls*.

    Fields are discovered by looking for descriptors that have both
    ``__get__`` and ``__set__`` (i.e. nanobind ``def_rw`` properties).
    Unset ``std::optional`` fields raise ``RuntimeError`` on access and
    are silently skipped.
    """
    def __repr__(self):
        fields = {}
        for name, val in cpp_cls.__dict__.items():
            if not name.startswith('_') and hasattr(val, '__get__') and hasattr(val, '__set__'):
                try:
                    fields[name] = getattr(self, name)
                except RuntimeError:
                    pass
        lines = [f"{cpp_cls.__name__}("]
        for k, v in fields.items():
            try:
                lines.append(f"  {k}={v!r},")
            except RuntimeError:
                lines.append(f"  {k}=<{type(v).__name__}>,")
        lines.append(")")
        return "\n".join(lines)
    return __repr__


def _numpy_setattr(cls_name):
    """Create a ``__setattr__`` that converts numpy types before delegating to C++.

    The returned method calls :func:`_to_builtin` on every value, then
    forwards to the C++ base-class setter.  On ``TypeError`` it re-raises
    with a message that includes the class name, field name, and value.
    """
    def __setattr__(self, name, value):
        value = _to_builtin(value)
        try:
            super(type(self), self).__setattr__(name, value)
        except TypeError as e:
            hint = " (must be a non-negative integer)" if isinstance(value, int) and value < 0 else ""
            raise TypeError(
                f"{cls_name}.{name} = {_fmt_val(value)} ({type(value).__name__}){hint}: {e}"
            ) from None
    return __setattr__


class InitialConditionSettings(_CppInitialConditionSettings):
    __doc__ = _CppInitialConditionSettings.__doc__
    __setattr__ = _numpy_setattr("InitialConditionSettings")


class OptimTargetSettings(_CppOptimTargetSettings):
    __doc__ = _CppOptimTargetSettings.__doc__
    __setattr__ = _numpy_setattr("OptimTargetSettings")


class ControlParameterizationSettings(_CppControlParameterizationSettings):
    __doc__ = _CppControlParameterizationSettings.__doc__
    __setattr__ = _numpy_setattr("ControlParameterizationSettings")


class ControlInitializationSettings(_CppControlInitializationSettings):
    __doc__ = _CppControlInitializationSettings.__doc__
    __setattr__ = _numpy_setattr("ControlInitializationSettings")


class Setup(_CppSetup):
    __doc__ = _CppSetup.__doc__
    __setattr__ = _numpy_setattr("Setup")

    __repr__ = _make_repr(_CppSetup)
