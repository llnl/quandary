"""Python subclass of Setup with __repr__ and improved error messages."""

from .._quandary_impl import Setup as _CppSetup


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
