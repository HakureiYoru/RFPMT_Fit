# Curvefit

This is the released version, the developing/debugging/preview version in: [Development Branch](https://github.com/HakureiYoru/curvefit/tree/differential_evolution)

## Simulation of electron motion in a fast changing field for new detector development

## Overview

This project involves the simulation of electron motion in rapidly changing fields, a process critical for the development of new detectors. At its core, the project is a data analysis and fitting application written in Python. It employs a set of parametric equations for fitting and analyzing calibration data and then applies the results to test data. The results are visualized on a map for easy interpretation and further analysis.

## Scripts

The project consists of three main Python scripts:

1. `UI.py`: This script provides a graphical user interface (GUI) that allows users to load data, set parameter boundaries, initiate a fitting process, view logs, and choose whether to autoscale graphs.

2. `Filter.py`: This script defines a set of parametric equations and uses these to fit the provided x and y data. The fitting process is performed in a separate thread to avoid blocking the GUI. The results of the fitting process are logged and saved to a file.

3. `function_analysis.py`: This script contains several functions used for further processing and analysis of the x and y data. Functions include Fourier transformations, keeping data within one period, and generating maps of the data.

## Workflow

1. **Calibration**: Load the calibration data (end0). Upon loading, a window displaying the pre-read parameter results will pop up.

2. **Fitting**: Click the "Fit" button to execute the fitting process. The result is a set of calibrated reference data.

3. **Testing**: After fitting, open a new "map" window. In this window, load the actual test data. The results of the test data will be displayed on various maps.

## Requirements

- Python 3.7 or later
- Libraries: numpy, scipy, Tkinter, json, logging

## Installation

1. Clone this repository to your local machine.
2. Install the necessary libraries by running `pip install -r requirements.txt` in your terminal.

## Usage

1. Run the `UI.py` script to launch the graphical user interface.
2. Load your data file, set parameter boundaries as needed.
3. Click "Fit" to start the fitting process.
4. After fitting, open the "map" window and load the test data.

## Algorithms

`function_analysis.py` contains several key algorithms and functions that warrant further explanation:

1. **Fourier Transform**: `xy_fft` function in `function_analysis.py` performs a Fourier Transform on the x and y data. Fourier Transform is a mathematical technique that transforms a function of time (a signal) into a function of frequency. This is used in the script to find the dominant frequency components of the provided data.

2. **Adjusting Frequency**: The `process_data` function adjusts the frequencies f1 and f2 based on the number of points in a cycle. This adjustment is crucial for subsequent analysis as it ensures that the frequency of the data is correctly represented.

3. **Data Segmentation**: The `keep_one_period` function is used to keep only one period of data. It calculates the Fourier Transform and finds the dominant frequency to determine the period of the data.

4. **Map Calculation**: The `calculate_map` function creates a map that describes the trajectory of the data. This map is created by taking into account the number of fitted points at each pixel, time information of each pixel, and the distance transform of the trajectory. It further computes a weight map from the distance transform and calculates the weights for each timestamp based on the time difference to the real time of the particle in this pixel.

These algorithms play a critical role in the simulation and analysis of electron motion in rapidly changing fields. Understanding these will aid in further development and application of this project.

## Usage Example

In the Calibration step, a typical data file that can be loaded might be a .dat file containing two columns of numbers, which represent the x and y values. Upon loading the data, the script will display the pre-read parameter results. You can then set the parameter boundaries as needed and initiate the fitting process by clicking the "Fit" button. 

After the fitting process, you can open a new "map" window. In this window, you can load the actual test data. The results of the test data will be displayed on various maps. These results can be used for further analysis and interpretation.

In all these steps, the scripts employ advanced algorithms such as Fourier Transform, adjusting frequency, data segmentation, and map calculation to ensure accurate and useful results.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

Please see the [LICENSE](./LICENSE) file for details.
