# pyemfield: Python Antenna Radiation Analysis and Optimization Library

`pyemfield` is a Python library designed for antenna radiation pattern analysis and gain optimization. The package leverages SciPy, NumPy, Matplotlib, and ANSYS HFSS (via `pyaedt`) to automate high-frequency structure simulations, process electromagnetic data, and visualize results.

## Key Features
- **Gain Calculation and Optimization**: Provides methods to calculate antenna gain with multiple optimization options, including exhaustive, heuristic, and SciPy-based optimizations.
- **Radiated Power Density and Cumulative Distribution**: Calculates and visualizes cumulative distribution functions (CDFs) for radiation characteristics.
- **Custom Radiation Pattern and Gain Distribution Plots**: Offers multiple plotting functions to visualize radiated power density, gain, and other related data.
- **Integration with HFSS**: Automates the export of far-field data from HFSS simulations for further analysis and optimization.

## Installation

You can install the package using `pip`:

```bash
pip install pyemfield
```

> **Note**: Ensure that `pyaedt` and ANSYS HFSS are installed to fully utilize all functionalities in this package.

## Dependencies
This package requires the following libraries:
- `numpy`
- `scipy`
- `matplotlib`
- `pyaedt`

## Quick Start

### 1. Basic Usage Example

The following example demonstrates how to use `pyemfield` to calculate and visualize antenna gain:

```python
from pyemfield import hfss_design, get_ffds, Beam

folder = r'D:\OneDrive - ANSYS, Inc\GitHub\pyemfield\tests\ffds'
ffds = get_ffds(folder)

ffds.keys()

x = {j:(1,0) for i, j in ffds.items()}
b1 = Beam(x)

b1.ffd_excitation
b2 = b1.optimize_gain(60, 60)
b2.plot_realized_gain_contour()
b2.ffd_excitation

hd.update_excitation(b2.ffd_excitation)
```

![2024-11-03_09-12-52](/assets/2024-11-03_09-12-52.png)

### 2. Create a Custom Plane and Plot Cumulative Distribution

```python
# Set up multiple beams and create a Plane
plane = create_plane("MyPlane", beam, [(0, 0), (90, 90)], fast=True)
plane.plot_rGain_cdf()
plane.plot_eirp_cdf()
```


