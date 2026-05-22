"""Regression test configuration for pytest."""


def pytest_addoption(parser):
    group = parser.getgroup("regression")
    group.addoption(
        "--exact",
        action="store_true",
        default=False,
        help="Use exact comparison for floating point numbers"
    )
