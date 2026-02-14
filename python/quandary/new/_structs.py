"""Python subclasses of nanobind-bound C++ structs with __repr__ and improved errors."""

from .._quandary_impl import (
    Setup as _CppSetup,
    InitialConditionSettings as _CppInitialConditionSettings,
    OptimTargetSettings as _CppOptimTargetSettings,
    ControlParameterizationSettings as _CppControlParameterizationSettings,
    ControlInitializationSettings as _CppControlInitializationSettings,
)


def _fmt_val(v):
    r = repr(v).replace('\n', ' ')
    return r[:77] + '...' if len(r) > 80 else r


def _make_repr(cpp_cls):
    """Create a __repr__ method that discovers fields from the C++ base class."""
    def __repr__(self):
        fields = {
            name: getattr(self, name)
            for name, val in cpp_cls.__dict__.items()
            if not name.startswith('_') and hasattr(val, '__get__') and hasattr(val, '__set__')
        }
        lines = [f"{cpp_cls.__name__}("]
        for k, v in fields.items():
            lines.append(f"  {k}={v!r},")
        lines.append(")")
        return "\n".join(lines)
    return __repr__


class Setup(_CppSetup):
    def __setattr__(self, name, value):
        try:
            super().__setattr__(name, value)
        except TypeError as e:
            hint = " (must be a non-negative integer)" if isinstance(value, int) and value < 0 else ""
            raise TypeError(
                f"Setup.{name} = {_fmt_val(value)} ({type(value).__name__}){hint}: {e}"
            ) from None

    __repr__ = _make_repr(_CppSetup)


class InitialConditionSettings(_CppInitialConditionSettings):
    __repr__ = _make_repr(_CppInitialConditionSettings)


class OptimTargetSettings(_CppOptimTargetSettings):
    __repr__ = _make_repr(_CppOptimTargetSettings)


class ControlParameterizationSettings(_CppControlParameterizationSettings):
    __repr__ = _make_repr(_CppControlParameterizationSettings)


class ControlInitializationSettings(_CppControlInitializationSettings):
    __repr__ = _make_repr(_CppControlInitializationSettings)
