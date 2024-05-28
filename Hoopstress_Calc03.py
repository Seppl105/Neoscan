# Calculating the hoopstresses with E(r) and nu(r) modelled as Fourier Series

import matplotlib.pyplot as plt
import numpy as np
import sympy as smp
from scipy.integrate import odeint
from scipy.integrate import solve_bvp


#####   Definierte Werte
anzahlRowPlots = 2
anzahlColumnPlots = 2

numberOfValues = 1000
#   Spannungsrandbedingungen
#s_z0 = 0    *  (10**6)  # [Pa] !!!Wert frei gewählt
#s_ra = 10   *  (10**6)  # [Pa] !!!Wert frei gewählt
#s_ri = 10   *  (10**6)  # [Pa] !!!Wert frei gewählt
s_z0 = 0
s_z = s_z0
ds_z = 0
s_r0 = 0
s_phi0 = 4.3264372825278825 * 10**8


#   Konstanten des Versuchsaufbau

r_i = 0.430 # [m] innerer radius
r_a = 0.646 # [m] äußerer radius
r = np.linspace(r_i,r_a, numberOfValues) # array mit diskreten Radien

j = 114.2   * (1000**2) # [A/m^2] Stromdichte
b_za = -2.5             # [T] magnetische Flussdichte außen
b_zi = 14               # [T] magnetische Flussdichte innen
b_0 = b_za - (b_za-b_zi)/(r_a-r_i) * r_a # [T] absolutes Glied der Geradengleichung für B(r)

def calcBFeld(r, r_a, r_i, b_za, b_zi, b_0):
    return ((b_za - b_zi)/(r_a - r_i))  *  r + b_0

#   Materialparameter
E = np.ones(numberOfValues) * 100 * 10**9 # E Modul
#E[:int(0.6* numberOfValues)] = 150* 10**9

#nu = 0.3             # [-] possion's ratio
#p = 1 - nu # [-] const.
nu = np.ones(numberOfValues) * 0.3 # Querkontraktionszahl
#nu[:int(0.3* numberOfValues)] = 0.25



# Fuktionswert an der Stelle r über eine Fourierreihe aus den Fourierkoeffizienten bestimmen
def inverseFourier(r, coefficients, frequencies):
    n = len(coefficients)
    result = np.real(coefficients[0]) / n # konstanter Anteil
    
    for k in range(1, n//2):  # aus symmetriegründen nur die erste Hälfte der Freuqenzen zu betrachten
        amplitude = 2 * np.abs(coefficients[k]) / n  # Amplitude berechnen
        phase = np.angle(coefficients[k])  # Phase berechnen
        result += amplitude * smp.cos(2 * 3.141592653589793 * frequencies[k] * r + phase)
        
    return result

# Fuktionswert an der Stelle r über eine Fourierreihe aus diskreten Inputwerten bestimmen
def getFourierSeries(r, inputFunction):
    fftInputFunction = np.fft.fft(inputFunction) # Koeffizienten der fft berechnet mit der Funktion fft aus dem fft package
    return inverseFourier(r, fftInputFunction, np.fft.fftfreq(len(inputFunction)))



######  Hauptprogramm

### Berechnen der differenzierbaren Funktionen E(r) und nu(r) und ihrer ersten Ableitungen nach r
# E(r)
plt.subplot(anzahlRowPlots, anzahlColumnPlots, 1)
plt.plot(r , E, label="vorgegebenes E(r)")
plt.xlabel("Radius in m")
plt.ylabel("E-Modul in N/m^2")

x = smp.Symbol("x", real=True)
fourierFunctionE = getFourierSeries(x, E) # A function dependant on the "Symbol" x
#print(FourierFunctionE)
fourierFunctionE_f = smp.lambdify(x, fourierFunctionE) # Convert FourierFunctionEr to a numerical function for plotting
#FourierFunctionEr = getFourierSeries(np.linspace(0, len(r), len(r)*5), Er)

plt.plot(np.linspace(r_i, r_a, numberOfValues*5), 
         fourierFunctionE_f(x=np.linspace(0, len(r), len(r) * 5)), label='mit Fourrierreihe genähert', linestyle='--')
plt.legend()

# nu(r)
plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlRowPlots + 1)
plt.plot(r , nu, label="vorgegebenes nu(r)")
plt.xlabel("Radius in m")
plt.ylabel("nu in 1")

fourierFunctionNu = getFourierSeries(x, nu) # A function dependant on the "Symbol" x
fourierFunctionNu_f = smp.lambdify(x, fourierFunctionNu) # Convert FourierFunctionEr to a numerical function for plotting
#FourierFunctionNu = getFourierSeries(np.linspace(0, len(r), len(r)*5), nu)
plt.plot(np.linspace(r_i, r_a, numberOfValues*5), 
         fourierFunctionNu_f(x=np.linspace(0, len(r), len(r) * 5)), label='mit Fourrierreihe genähert', linestyle='--')
plt.legend()

# dE(r)/dr
dfourierFunctionE = smp.diff(fourierFunctionE, x)
dfourierFunctionE_f = smp.lambdify(x, dfourierFunctionE)

# dnu(r)/dr
dfourierFunctionNu = smp.diff(fourierFunctionNu, x)
dfourierFunctionNu_f = smp.lambdify(x, dfourierFunctionNu)


### Bestimmen der Spannungen in radiale Richtungen s_r(r) und in Umfangsrichtung s_phi(r)

# calculating the radial force, the numerical pendant, the first derivative and the numerical pendant
radialForce = calcBFeld(x, r_a, r_i, b_za, b_zi, b_0) * j
radialForce_f = smp.lambdify(x, radialForce)
dradialForce = smp.diff(radialForce, x)
dradialForce_f = smp.lambdify(x, dradialForce)

# lösen des DGL Systems 1. Ordnung bestehend aus (b*) für d(sigma_r)/dr und (a*) für d(sigma_phi)/dr

def dSdr(r, S):
    s_r, s_phi = S
    return [    1/r * s_phi   - calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   - 1/r * s_r,
                1/r * s_phi   - calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   - 1/r * s_r   + fourierFunctionNu_f(r) * ds_z   + dfourierFunctionNu_f(r) * (s_r + s_z)   + dfourierFunctionE_f(r) * 1/fourierFunctionE_f(r) * (-fourierFunctionNu_f(r) * s_r  + s_phi  + fourierFunctionNu_f(r) * s_z)   - (1 + fourierFunctionNu_f(r)) * radialForce_f(r)]

## solving initial value problem with odeint
# defining initial conditions for s_r(r_i) and s_phi(r_i)   (initial conditions are the values for r[0])
S_0 = (s_r0, s_phi0)
solOdeint = odeint(dSdr, y0=S_0, t=r, tfirst=True)


s_rsolOdeint = solOdeint.T[0]
s_phisolOdeint = solOdeint.T[1]


plt.subplot(anzahlRowPlots, anzahlColumnPlots, 2)
plt.plot(r, s_rsolOdeint, label="s_r with Odeint")
plt.xlabel("Radius in m")
plt.ylabel("s_r in N/m^2")
plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlRowPlots + 2)
plt.plot(r, s_phisolOdeint, label="s_phi with Odeint")
plt.xlabel("Radius in m")
plt.ylabel("s_phi in N/m^2")

# ## solving bvp with solve_bvp
# def bc(ya, yb):
#     return np.array([ya[0], ya[-1]])

# Sguess = np.zeros((len(r), len(r)), dtype=float)
# solBvp = solve_bvp(dSdr, bc, r, Sguess)

# s_rsolBvp = solBvp.T[0]
# s_phisolBvp = solBvp.T[1]


# plt.subplot(anzahlRowPlots, anzahlColumnPlots, 2)
# plt.plot(r, s_rsolBvp, label="s_r with bvp")
# plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlRowPlots + 2)
# plt.plot(r, s_phisolBvp, label="s_phi with bvp")

# # lösen der Gleichung (e) also eihner DGL 2. Ordnung durch überführen in ein DGL System mit DGLs 1. Ordnung

# # Vektor containing s_r and ds_r/dr = h
# def dSdr(r, S):
#     s_r, h = S
#     return [    h,
#                 -1/r * 3 * h    - dradialForce_f(r)    - 1/r * (+fourierFunctionNu(r) * ds_z + dfourierFunctionNu_f(r) * (s_r + s_z) + dfourierFunctionE_f(r) * 1/fourierFunctionE(r) * (-fourierFunctionNu(r) * s_r + r * h + r*radialForce_f(r) + s_r - fourierFunctionNu(r) * s_z) - (2 + fourierFunctionNu_f(r)) * radialForce_f(r))] 
#

# # defining initial conditions for s_r(r_i) and ds_r(r_i)   (initial conditions are the value for r[0])
plt.show()