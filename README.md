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
from pyemfield import Ffd, Beam, Plane, hfss_design

# Initialize HFSS design and export FFD files
hd = hfss_design()
folder = "path_to_ffd_folder"
hd.export_ffds(folder, 'Setup1 : Sweep', '30.0GHz')

# Load FFD files
ffds = get_ffds(folder)

# Create a Beam object and optimize gain
beam = Beam({ffd_obj: (1, 0) for ffd_obj in ffds.values()})
optimized_beam = beam.optimize_gain(40, 60)

# Visualize the gain contour plot
optimized_beam.plot_realized_gain_contour()
```

### 2. Create a Custom Plane and Plot Cumulative Distribution

```python
# Set up multiple beams and create a Plane
plane = create_plane("MyPlane", beam, [(0, 0), (90, 90)], fast=True)
plane.plot_rGain_cdf()
plane.plot_eirp_cdf()
```

## Contact
If you have any questions or feedback, feel free to reach out:
- **Email**: mingchih.lin8@gmail.com
- **GitHub**: [pyemfield Repository](https://github.com/linmingchih/pyemfield)

## License
This project is licensed under the MIT License. For more details, see the [LICENSE](https://opensource.org/licenses/MIT) file.

