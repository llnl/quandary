"""Python subclasses of C++ structs with numpy-to-builtin conversion."""

import numpy as np

from .._quandary_impl import Setup as _CppSetup
from .._quandary_impl import (
    InitialConditionSettings as _CppInitialConditionSettings,
    OptimTargetSettings as _CppOptimTargetSettings,
    ControlParameterizationSettings as _CppControlParameterizationSettings,
    ControlInitializationSettings as _CppControlInitializationSettings,
)


def _to_builtin(value):
    """Convert numpy arrays/scalars to Python builtins for nanobind."""
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
    r = repr(v).replace('\n', ' ')
    return r[:77] + '...' if len(r) > 80 else r


def _make_repr(cpp_cls):
    """Create a __repr__ method that discovers fields from the C++ base class."""
    def __repr__(self):
        fields = {}
        for name, val in cpp_cls.__dict__.items():
            if not name.startswith('_') and hasattr(val, '__get__') and hasattr(val, '__set__'):
                try:
                    fields[name] = getattr(self, name)
                except RuntimeError:
                    # Skip unset optional fields (bad_optional_access)
                    pass
        lines = [f"{cpp_cls.__name__}("]
        for k, v in fields.items():
            try:
                lines.append(f"  {k}={v!r},")
            except RuntimeError:
                # Handle nested objects with unset optional fields
                lines.append(f"  {k}=<{type(v).__name__}>,")
        lines.append(")")
        return "\n".join(lines)
    return __repr__


def _numpy_setattr(cls_name):
    """Create a __setattr__ that converts numpy types before delegating to C++."""
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
    __setattr__ = _numpy_setattr("InitialConditionSettings")


class OptimTargetSettings(_CppOptimTargetSettings):
    __setattr__ = _numpy_setattr("OptimTargetSettings")


class ControlParameterizationSettings(_CppControlParameterizationSettings):
    __setattr__ = _numpy_setattr("ControlParameterizationSettings")


class ControlInitializationSettings(_CppControlInitializationSettings):
    __setattr__ = _numpy_setattr("ControlInitializationSettings")


class Setup(_CppSetup):
    """Python subclass of the C++ Setup struct.

    Extends the nanobind-bound C++ base class with:
    - A readable ``__repr__`` that shows all field values.
    - Improved ``TypeError`` messages with a hint when a negative integer
      is assigned to an unsigned field.
    - Automatic conversion of numpy types to Python builtins for nanobind.
    """
    __setattr__ = _numpy_setattr("Setup")

    __repr__ = _make_repr(_CppSetup)
