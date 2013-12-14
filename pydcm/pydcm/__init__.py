# Python module pydcm
# Create straightforward applications for simulations using coolfluid3
#
# Space discretisation uses the Discontinuous Collocation Methods
# - Spectral Difference Method
#
# A wide range of Time discretisation exists:
# - Explicit Runge Kutta methods
# - Semi-implicit LUSGS methods (BDF1, BDF2)
#
# Implemented equations are:
# - Linear Advection Diffusion
# - Euler
# - Navier-Stokes 


# Bring specific nested scopes up
from advection_diffusion   import AdvectionDiffusion
from euler                 import Euler
from navierstokes          import NavierStokes
from lineuler              import LinearizedEuler
from case                  import save, restart
from util                  import log
