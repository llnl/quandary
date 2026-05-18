# Performance Regression Test Documentation

These tests are run with `pytest` and the pytest-benchmark plugin.

To run tests from base directory:
```
pytest tests/performance/
```
Or from current directory:
```
pytest
```

## Useful options:

- `pytest -v` will print the names of the tests being run

- `pytest -s` will print output to screen

- `pytest -k "myConfig_4"` runs only the test with name "myConfig" and with 4 processes

See `pytest --help` for more options.

## Local performance measurement
To save performance data from a run do:
```
pytest -s -k "test_of_interest_4" --benchmark-autosave
```
This saves results to the `.benchmarks` directory.

To compare a current run to the previous run do:
```
pytest -s -k "test_of_interest_4" --benchmark-compare
```

For more options, see pytest-benchmark [documentation](https://pytest-benchmark.readthedocs.io/en/stable/comparing.html).

## How to add a test

1. Create a config file in the `configs` directory, e.g., `tests/performance/configs/newSimulation.toml`
2. Add new entry to test_cases.json:
    - The `simulation_name` should be the new simulation name, e.g. `newSimulation`.
    - The `number_of_processes` is an array of integers. For each integer `i`, a simulation will be run with `mpirun -n ${i}`.
    - The `repetitions` is the number of times to run this test and to average the timings over.
