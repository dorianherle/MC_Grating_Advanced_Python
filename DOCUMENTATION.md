# MC Grating Python Documentation

This project is a Python wrapper for the MC Grating software. It provides a Pythonic interface to the software's functionality, allowing users to create geometries, run simulations, and visualize the results all within Python.

## File Structure

The project consists of several Python files, each containing different classes and functions:

1. `run_simulation.py`: Contains functions for running a simulation and processing the results.
2. `visualizations.py`: Contains classes for visualizing the geometry and the results of a simulation.
3. `mc_grating_python.py`: Contains the main class of the program, `mc_grating_python`, which provides methods for creating geometries, running simulations, and visualizing the results.

## Usage

To use this project, you first need to create a geometry using the `geo` method of the `mc_grating_python` class. This method takes the periods in the x and y directions as arguments and returns a `geometry3D` object.

Next, you can visualize the geometry using the `show` method of the `mc_grating_python` class. This method takes a `geometry3D` object and an optional color argument, and returns a `visualization` object.

To run a simulation, use the `sim` method of the `mc_grating_python` class. This method takes a `geometry3D` object and several optional arguments for the simulation parameters, and returns an `analyze` object.

Finally, you can plot the results of the simulation using the `plot` method of the `mc_grating_python` class. This method takes a result object, and optional geometry and color arguments, and returns a `plot_res` object.

## Example

Here is an example of how to use this project:

```python
from mc_grating_python_ver1.mc_grating_python import mc_grating_python

# Create a geometry
geometry = mc_grating_python.geo(1, 1)

# Visualize the geometry
visualization = mc_grating_python.show(geometry)

# Run a simulation
analyze = mc_grating_python.sim(geometry)

# Plot the results
plot = mc_grating_python.plot(analyze)
```
