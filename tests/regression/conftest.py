"""Regression test configuration for pytest."""


def pytest_addoption(parser):
    group = parser.getgroup("regression")
    group.addoption(
        "--exact",
        action="store_true",
        default=False,
        help="Use exact comparison for floating point numbers"
    )
    group.addoption(
        "--config-format",
        action="store",
        choices=["toml", "cfg"],
        default="toml",
        help="Configuration file format to use (default: toml)"
    )
