## Stochastic Simulations of Neuronal Avalanches on Wilson-Cowan Based Population Models

This repo is an implementation of widely used stochastic simulation algorithms (SSA, T-Leaping, AdT-Leaping) for Wilson-Cowan based population models for subsequent analyses of neuronal avalanches and power law correlations. For preceding research on stochastic Wilson-Cowan models and neuronal avalanches refer to the below articles respectively.
<br>
```
Adrián Ponce-Alvarez, Adrien Jouary, Martin Privat, Gustavo Deco, Germán Sumbre,
Whole-Brain Neuronal Activity Displays Crackling Noise Dynamics,
Neuron,
Volume 100, Issue 6,
2018,
Pages 1446-1459.e6,
ISSN 0896-6273,
https://doi.org/10.1016/j.neuron.2018.10.045.
```

```
Adrián Ponce-Alvarez, Adrien Jouary, Martin Privat, Gustavo Deco, Germán Sumbre,
Whole-Brain Neuronal Activity Displays Crackling Noise Dynamics,
Neuron,
Volume 100, Issue 6,
2018,
Pages 1446-1459.e6,
ISSN 0896-6273,
https://doi.org/10.1016/j.neuron.2018.10.045.
```
<br>

Wilson-Cowan models are defined by a population of neurons that are either in quiescent or active state with varying transition dynamics between these 2 states. On a stochastic model, these dynamics are interpreted probabilistically which coincides with the "Master Equation" in chemical statistics. This probabilistic explanation allows borrowing ideas from chemical kinematics and harnessing stochastic simulators for simulating continuous and deterministic dynamics through a discrete and stochastic simulator. Here, we address these simulations through Gillespie's Algorithm and Tau Leaping Algorithm in which the latter is event based and exact whereas the former iterates over a fixed time scheduler and is therefore approximate.
<br>
<br>
Stochastic Wilson-Cowan models, which show avalanche-like temporal and spatial spikes, can be analyzed with varying degrees of connectivity as candidates models for population dynamics in analyzing neuronal avalanches. Here, of course, the definition of spike bins dictates the accuracy and physical realness of the model. In this implementation we compare 3 different methods, namely, non-spatial spike counting, non-spatial Calcium events and spatio-temporal Calcium events.

## System Requirements
Simulations are all tested and run on MATLAB R2025b with NVIDIA RTX4070 GPU and intel i-14650HX CPU. For this setup, simulation times change between 30 mins to 5 hours depending on the available data.

## Running the Simulations
To run the simulations just execute the main program run_EImodel.m. Input and simulator types are hard coded, if necessary they can indeed be changed inside the main program by referring to the related sections.

### Inputs
For completeness we consider 3 different input types: constant(0.001 mag), noisy stepwise jumps around constant mean(0.001 mean, 0.0005 var) and sinusoidal (0.001 amp, 0.1 freq). It is possible to switch between these three input paradigms by changing the input mode, refer to:

```
% Inputs:
inputmethod = 0;
if inputmethod == 0
    I = 0.001*ones(N,t_max);
elseif inputmethod == 1
    mean = 0.001;
    var = const*0.5;
    I = var*(rand([N,t_max+5])*2-1) + mean;
elseif inputmethod == 2
    amp = 0.001;
    freq = 0.1;
    I = amp*sin(freq*t_min:0.01:t_max);
end
```

### Simulators
We use naive versions of Gillespie's Algorithm and Tau Leaping Algorithm to evaluate the temporal evaluation of the model. It is possible to switch between these two models by changing the operation mode, refer to:

```
% Simulator:
simmethod = 0;
if simmethod == 0
    simulator = Gillespie_EImodel;
    args = {};
elseif simmethod == 1
    simulator = TauLeaping_EImodel;
    args = {mode:0};
elseif simmethod == 2
    simulator = TauLeaping_Adaptive_EImodel;
    argS = {mode:0, base:0};
end
```

## Citing

If you find this work useful, please consider citing it.

```
@article{forgione2023from,
  author={Forgione, Marco and Pura, Filippo and Piga, Dario},
  journal={IEEE Control Systems Letters}, 
  title={From System Models to Class Models:
   An In-Context Learning Paradigm}, 
  year={2023},
  volume={7},
  number={},
  pages={3513-3518},
  doi={10.1109/LCSYS.2023.3335036}
}
```
## License

This repository is released under the MIT license. See [LICENSE](LICENSE) for additional details.