
## Discontinuous Collocation Methods ##

This repository contains two plugins that provide extra functionality to the coolfluid 3 project.
- DCM:  Discontinuous Collocation Methods
- SDM:  Spectral Difference Method

**Plugin DCM**

Discontinuous Collocation Methods acts as a library of building blocks
for a family of discontinuous collocation methods such as the 
spectral difference method, discontinuous galerkin method, ...

Provided in this plugin are e.g. libraries
- core: Metrics (interpolation, etc), element connectivities, ...
- equations: Advection Diffusion, Euler, Linearized Euler, Navier-Stokes
- solvers: classic Runge-Kutta, Low-storage Runge-Kutta

**Plugin SDM**

Spectral Difference Method for solving Partial Differential Equations

Provided in this plugin are, using building blocks from the DCM plugin:
- Computation of convective and diffusive terms using the Spectral Difference
  method. Diffusive terms are computed using the Second approach of Bassi-Rebay (BR2)
- Low-storage explicit Runge-Kutta schemes, specifically optimised for the SDM

###Installation: ###

  + Create a plugins directory if not existing, and clone this repository inside:

```
mkdir -p $PLUGIN_DIR
cd $PLUGIN_DIR
git clone https://github.com/coolfluid/sdm.git $PLUGIN_DIR/dcm_plugins
```

  + Rerun cmake in the coolfluid3 build directory:

```
cd $CF3_BUILD_DIR
cmake .  -DCF3_PLUGIN_DIRS=$PLUGIN_DIR
```
