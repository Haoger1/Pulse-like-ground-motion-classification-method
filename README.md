# Pulse-like-ground-motion-classification-method

This is the algorithm for the paper entitled "A new method for the identification and classification of pulse-like ground motions". This algorithm can identify and classify different pulse-like ground motion according ot identified pulse number, e,g., single pulse-like ground motion, multiple pulse-like ground motion, etc. The method is keeping improved and debugged even after article acceptance. Please read the paper for more details.

## Disclaimer
This algorithm is distributed under the GNU General Public License v3.0. Please refer to the LICENSE file for details.

## Algorithm

This repository provides the MATLAB implementation of the proposed pulse-like ground-motion identification and classification method.

Before using this classification algorithm, users should first apply the **Shock Waveform Decomposition (SWD)** method to the original velocity time series. The SWD procedure decomposes each ground motion record into several shock-waveform (SW) components and generates the variable `SWresults`, which is used as the input of the present classification algorithm. The SWD method has been made publicly available and can be accessed at [Shock Waveform Decomposition Method](https://github.com/galois-yan/Shock_Waveform_Decomposition_Method).

### Files

This repository contains two main MATLAB files:

```text
classifyPulselikeGM.m
classifyPulselikeGM_example_run.m
```

- `classifyPulselikeGM.m` is the main classification function.
- `classifyPulselikeGM_example_run.m` is an example script showing how to call the function.

### Input data

The input variable should be named `SWresults`.

`SWresults` should be a cell array. The first column of each cell stores the decomposed results obtained from SWD. Each data matrix in `SWresults{i,1}` should follow the structure below:

```text
Column 1   : Time column of the original time series
Column 2   : Velocity column of the original time series
Column 3   : Energy proportion of each decomposed SW component relative to the original signal
Column 4   : First decomposed SW component
Column 5   : Second decomposed SW component
Column 6   : Third decomposed SW component
...
Column 3+j : The j-th decomposed SW component, where j = 1, 2, ..., n
```

For Column 3, each valid row corresponds to one decomposed SW component. Therefore, if a record is decomposed into `n` SW components, the first `n` valid rows in Column 3 store the corresponding energy proportions of these components.

### Basic usage

Put `classifyPulselikeGM.m`, `classifyPulselikeGM_example_run.m`, and your `SWresults` data file in the same folder. Then open MATLAB and run:

```matlab
clear; clc;

load('SWresults.mat');  % This file should contain the variable SWresults.

results = classifyPulselikeGM(SWresults, true);

disp(struct2table(results));
```

The second input argument controls whether the classification information is displayed in the MATLAB Command Window:

```matlab
results = classifyPulselikeGM(SWresults, true);   % Display classification messages
results = classifyPulselikeGM(SWresults, false);  % Do not display classification messages
```

### Example script

Users can also directly run the example script:

```matlab
classifyPulselikeGM_example_run
```

This script loads the example `SWresults_sample` data, calls the main classification function, and displays the classification results.

### Output

The function returns a structure array named `results`. Each element corresponds to one ground-motion record.

The main fields include:

```text
recordIndex              Index of the ground-motion record
class                    Final classification result: SPGM, MPGM, or NPGM
subclass                 More detailed subclassification, such as DPGM
selectedComponents       Indices of selected decomposed SW components
firstEnergyProportion    Energy proportion of the first decomposed SW component
message                  Classification message displayed in MATLAB
```

The classification categories are:

```text
SPGM : Single pulse-like ground motion
DPGM : Double pulse-like ground motion
MPGM : Multiple pulse-like ground motion with more than two pulses (denoted here only since MPGM includes DPGM in original paper)
NPGM : Non-pulse-like ground motion
```

In this implementation, DPGM is treated as a special case of MPGM when two decomposed SW components are selected as pulse components.

### Requirements

This code requires MATLAB and the Signal Processing Toolbox, because the function `findpeaks` is used to detect local maxima of decomposed SW components.

### Notes

The classification algorithm assumes that the SWD results have already been obtained. This repository does not perform the SWD decomposition itself; it only identifies and classifies pulse-like ground motions based on the decomposed SW components stored in `SWresults`.

If your data use `NaN` instead of `0` to pad unused rows in Column 3, please modify the corresponding line in `classifyPulselikeGM.m` accordingly.

## Citation

If you use this code or the classification method in your research, please cite the following paper:

```text
Wang Junhao, Li Q.M. 2026. A new method for the identification and classification of pulse-like ground motions. Major revision
```

## Contact

For questions about algorithm and classification method, please contact: junhao.wang-3@postgrad.manchester.ac.uk
