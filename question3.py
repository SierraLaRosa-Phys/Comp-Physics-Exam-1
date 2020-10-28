from numpy import linspace, max
from toolbox import Jacobi_bc, bc_diff_eq
import matplotlib.pyplot as mp

def equation(n=5000, s=Jacobi_bc):
    p = lambda x: 1.0
    q = lambda x: 0.0

    coeff = [q,q,p]
    a = -1
    b = 1
    dx = (b-a)/(n+1)

    bc = {0:0.0, n:0.0}
    def rhs(x):
        if x == 0:
            return -1.0
        else:
            return 0.0

    x = linspace(a,b,n+1)
    y = bc_diff_eq(coeff, rhs, a, b, bc, n, s)

    return x, y, dx

position, voltage, delta = equation()

epsilon_0 = 8.85e-12
q_electron = -1.6e-19
delta_cubed = delta**3
phi_scale = q_electron/(epsilon_0*delta_cubed)

print(phi_scale*max(voltage))


position *= 2 #nm
voltage *= phi_scale*max(voltage)
mp.title(r"$\frac{d^2 \phi}{d x^2} = \frac{- \rho\left(x\right)}{\epsilon_0}$")
mp.ylabel(r"Potential $[V]$")
mp.xlabel(r"Position $[nm]$")
mp.plot(position,voltage, label="Potential")
mp.axvline(x=0, linestyle="--", label=r"Position of $q_e$")
mp.legend()
mp.show()
