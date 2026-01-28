
# Interfacing to Python environment
The python interface eases the use of the C++ Quandary code, and adds additional functionality: 

* Automatic estimation of the required time-step size based on eigenvalue decomposition of the system Hamiltonian
* Automatic computation of carrier frequencies based on system transition frequencies (Hamiltonian eigenvalue differences)
* Simulate and optimize with **custom system and control Hamiltonian operators**, other than the default model for superconducting quantum devices defined in [User Guide](user_guide.md). 

All interface functions are defined in `./quandary.py`, which defines the `Quandary` dataclass that gathers all configuration options and sets default values. Default values are overwritten by user input either through the constructor call through `Quandary(<membervar>=<value>)` directly, or by accessing the member variables after construction and calling `update()` afterwards. A good place to get started is to look at the example `example_swap02.py`, or jump start with the [Jupyter Notebook Tutorial](QuandaryWithPython_HowTo.ipynb)

## Parameter Naming

For historical reasons, some parameter names differ between the Python interface and TOML configuration. For a complete cross-reference of all parameters across interfaces, see the [Parameter Reference](../parameter_reference.md).

Under the hood, the python interface dumps all configuration options to file and evokes (multiple) subprocesses that execute the C++ code on this configuration file. The C++ output files are loaded back into the python interpreter.
It is therefore recommended to utilize the python interface only for smaller system sizes (fewer qubits), and switch to operate the C++ code directly when larger systems are considered (e.g. when parallel linear algebra is required.)

### Custom Hamiltonian models 
To enable custom Hamiltonians, pass the option `standardmodel=False` to the Quandary object and provide the complex-valued system Hamiltonian $H_d$ with `Hsys=<yourSystemHamiltonian>` as well as the real and imaginary parts of the custom control Hamiltonians per oscillator with `Hc_real=[<HcReal oscillator1, HcReal oscillator2, ...]` (will be multiplied by controls $p^k(t)$) `Hc_imag=[<HcImag oscillator1, HcImag oscillator2, ...]` (will be multiplied by controls $iq^k(t)$).

  * The units of the system Hamiltonian should be angular frequency (multiply $2\pi$), whereas the control Hamiltonian operators should be 'unit-free', since those units come in through the multiplied control pulses $p$ and $q$.
  * The control Hamiltonian operators are optional, but the system Hamiltonian is always required if `standardmodel=False`. 
  * The matrix-free solver can not be used when custom Hamiltonians are provided. The code will therefore be slower.

# Summary of all python interface options
Here is a list of all options available to the python interface, including their default values. Use `help(Quandary)` to get additional information on available interface functions, such as `quandary.simulate(...)`, `quandary.optimize(...)`.

```python
--8<-- "quandary.py:13:100"
```
