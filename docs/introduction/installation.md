# Installation

optmat requires Python 3.12 or later, NumPy 2.3 and SciPy 1.18. Install it from PyPI:

```sh
python -m pip install "optmat[gas]"
```

or from the [source repository](https://github.com/LupoLab/optmat), as below. The
package includes its
curated records and stored models; using them does not require the curators' paper
library or a new fit.

## Install with uv

With [uv](https://docs.astral.sh/uv/) installed:

```sh
git clone https://github.com/LupoLab/optmat.git
cd optmat
uv sync --locked --no-dev --extra gas
uv run --no-sync python
```

The last command opens Python in the environment just installed. Use
`uv run --no-sync python your_script.py` to run a script in that environment.

The `gas` extra supplies CoolProp for real-gas densities. It is needed
for the gas-index examples in the {doc}`/quickstart`, which use real-gas density by
default. For work that does not need this extra, omit `--extra gas`.

## Install with pip

Alternatively, clone the repository and install it in a virtual environment:

```sh
git clone https://github.com/LupoLab/optmat.git
cd optmat
python3 -m venv .venv
. .venv/bin/activate
python -m pip install ".[gas]"
python
```

These activation commands are for a POSIX shell. In Windows PowerShell, create the
environment with `py -3 -m venv .venv` and activate it with
`.venv\Scripts\Activate.ps1`. The remaining pip commands are the same. Use
`python -m pip install .` if real-gas densities are not needed.

## Check the installation

In Python, read the installed version:

```python
import optmat

optmat.__version__
```

Continue with the {doc}`/quickstart` for gas and solid indices, nonlinear
coefficients, Raman response, and ionisation rates. Record the package version and
the {doc}`identity of the models </identity>` when saving calculation results.

To run tests, build the manual, or change the package, use the development
environment described in the {doc}`/developers/index`.
