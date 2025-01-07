# Hydrograd
### ***Computational Hydrodynamics with Automatic Differentiation and Scientific Machine Learning***

A Julia-based computational hydrodynamics package with automatic differentiation (AD) capabilities. The package is designed to solve the partial differential equations (PDEs) of hydrodynamics, including the shallow water equations (SWE) and the Navier-Stokes equations (NS). Currently, the package is under development and only supports the SWE.

Hydrograd heavily relies on [SciML.jl](https://github.com/SciML/SciML.jl), which is a collection of tools for scientific machine learning and differential equation solving. 

## Overview
Hydrograd provides ***differentiable models*** for solving the hydrodynamics governing equations with finite volume methods. It supports:
- Forward simulation of flow dynamics
- Parameter inversion 
- Sensitivity analysis through automatic differentiation

For parameter inversion and sensitivity analysis, the package uses the automatic differentiation capabilities of [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and [Zygote](https://github.com/FluxML/Zygote.jl). The parameters currently supported are Manning's coefficient, bed elevation, and inlet discharge, which cover many common cases in the field of hydrodynamics, such as river flow, flood wave, and dam break. More parameters will be supported in the future. Users can also define their own parameters and use them in the parameter inversion and sensitivity analysis.

## Features
- The models are differentiable with respect to the parameters and the initial conditions. 
- Differentiable models can bridge physics-based and data-driven models, for example, the use of [universal differential equations (UDEs)](https://arxiv.org/abs/2001.04385) to bridge the gap between the shallow water equations and the neural networks.
- The models are designed to be compatible with the [SciML.jl](https://github.com/SciML/SciML.jl) ecosystem, which includes a wide range of tools for scientific machine learning and differential equation solving.
- Automatic differentiation using [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and [Zygote](https://github.com/FluxML/Zygote.jl)

## Other notable features
- Finite volume discretization
- Gudonov-type Riemann solver
- Support for various boundary conditions

- For 2D cases, the package directly supports the format of the [SRH-2D](https://www.usbr.gov/tsc/techreferences/computer%20software/models/srh2d/index.html) solver (USBR's Sedimentation and River Hydraulics - Two-Dimensional model; GUI with [SMS](https://aquaveo.com/software/sms/introduction)).

## Prerequisites
- Julia 1.11 or later
- Python 3.10 or later (for SRH-2D support)
- [pyHMT2D](https://github.com/psu-efd/pyHMT2D) package (for SRH-2D support)

## Installation

1. Clone the repository:

```bash
git clone https://github.com/psu-efd/Hydrograd.jl.git
cd Hydrograd.jl
```

2. Activate and instantiate the Julia environment:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This will automatically install all required Julia dependencies specified in `Project.toml`.

3. Run an example: See [Running Instructions](examples/SWE_2D/README.md) for detailed usage.

## License
MIT License

## Contact
Xiaofeng Liu, Ph.D., P.E.

Department of Civil and Environmental Engineering

Penn State University

Web: [https://water.engr.psu.edu/liu/](https://water.engr.psu.edu/liu/)

Email: [xiaofengliu19@gmail.com](mailto:xiaofengliu19@gmail.com)

For questions and support, please open an issue in the repository.
