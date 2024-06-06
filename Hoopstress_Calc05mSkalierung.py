import matplotlib.pyplot as plt
import numpy as np
import sympy as smp
from scipy.integrate import odeint
from scipy.integrate import solve_bvp
from functions.Hoopstress_Calc02Copy import calcStresses
from functions.calcMaterialFunction import calcMaterialTanhCoefficients
from functions.calcMaterialFunction import calcTanhValue
from datetime import datetime # Nur für die Bennenung der Grafiken

######## Bemerkungen
# Einheiten in SI Basiseinheiten mit mSkalierung kann die Längeneinheit skaliert werden
# 

#####   Definierte Werte
mSkalierung = 10**(3)        # Skalierungsfaktor für die Einheit Meter gleich 1, für mm gleich 10*(3), etc. 
windingDivisor = 60          # Es wird für 600/windingDivisor Windungen gerechnet
numberOfValues = int(100000 / windingDivisor) # Anzahl der Werte des r Vektors

#   Subplots
anzahlRowPlots = 2
anzahlColumnPlots = 4

#   Spannungsrandbedingungen

#s_z0 = 0    *  (10**6)  # [Pa] !!!Wert frei gewählt
s_ra = 0   *  (10**6)  # [Pa] !!!Wert frei gewählt
s_ri = 0   *  (10**6)  # [Pa] !!!Wert frei gewählt
s_z0 = 0
s_z = s_z0
ds_z = 0
#s_r0 = 0
#s_phi0 = 4.3264372825278825 * 10**8 * mSkalierung**(-1)     # N/m^2 E-3 = kg m/(s^2 m^2) E-3 = kg/(s^2 m) E-3 = kg/(s^2 mm)


#   Materialparameter

t = 0.36 * 10**(-3) * mSkalierung     # [m] Dicke einer Wicklung (Conductor, Cowinding, Insulation)
t_con = 0.12 * 10**(-3) * mSkalierung # [m] Dicke des Bandleiters
t_cow = 0.23 * 10**(-3) * mSkalierung # [m] Dicke des Cowindings
t_ins = 0.01 * 10**(-3) * mSkalierung # [m] Dicke der Isolation
materialWidths = [t_con, t_cow, t_ins]

E_con = 280 * 10**9 * mSkalierung**(-1) # E-Modul Conductor
E_cow = 300 * 10**9 * mSkalierung**(-1) # E-Modul Cowinding
E_ins = 200 * 10**9 * mSkalierung**(-1) # E-Modul Insulation
materialEs = [500, 450, 400]

ny_con = 0.35 * 10**9 * mSkalierung**(-1) # Possion's Ratio Conductor
ny_cow = 0.3 # Possion's Ratio Cowinding
ny_ins = 0.4 # Possion's Ratio Insulation
materialNys = [0.35, 0.3, 0.4]

#   Konstanten des Versuchsaufbau

r_i = 0.430 * mSkalierung               # [m] innerer radius
r_a = 0.646 * mSkalierung               # [m] äußerer radius
r_a = (r_i + t * ( ((r_a - r_i) / t) / windingDivisor))

r = np.linspace(r_i,r_a, numberOfValues) # array mit diskreten Radien

j = 114.2 * 10**6 * mSkalierung**(-2)       # [A/m^2] Stromdichte
b_za = -2.5                             # [T] magnetische Flussdichte außen
b_zi = 14                               # [T] magnetische Flussdichte innen
b_0 = b_za - (b_za-b_zi)/(r_a-r_i) * r_a # [T] absolutes Glied der Geradengleichung für B(r)

# Funktion zur Berechnung des B-Felds an der Stelle r
def calcBFeld(r, r_a, r_i, b_za, b_zi, b_0):
    return ((b_za - b_zi)/(r_a - r_i))  *  r + b_0

coefficientsE = calcMaterialTanhCoefficients(r_i , r_a, t, materialWidths, materialEs, slope=1000, scale=720)
E = calcTanhValue(r, coefficientsE, materialEs)
print(E)
#E = np.ones(numberOfValues) * 100 * 10**9 * mSkalierung**(-1) # E Modul in  N/m^2 E-3 = kg m/(s^2 m^2) E-3 = kg/(s^2 m) E-3 = kg/(s^2 mm)
#E[:int(0.2* numberOfValues)] = 150 * 10**9 * mSkalierung**(-1)
#E[len(E) - int(0.2* numberOfValues):] = 150 * 10**9 * mSkalierung**(-1)


#ny = 0.3             # [-] possion's ratio
#p = 1 - ny # [-] const.

#ny = np.ones(numberOfValues) * 0.3 # Querkontraktionszahl
coefficientsNy = calcMaterialTanhCoefficients(r_i, r_a, t, materialWidths, materialNys, slope=1000, scale=720)
ny = calcTanhValue(r, coefficientsNy, materialNys)
#ny[len(E) - int(0.2* numberOfValues):] = 0.25
#ny[:int(0.2* numberOfValues)] = 0.25



####################################################################################################################
##### Funktionsdefinition


# Fuktionswert an der Stelle r über eine Fourierreihe aus den Fourierkoeffizienten bestimmen
def inverseFourier(r, coefficients, frequencies):
    n = len(coefficients)
    result = np.real(coefficients[0]) / n # konstanter Anteil
    
    for k in range(1, n//2):  # aus symmetriegründen nur die erste Hälfte der Freuqenzen zu betrachten
        amplitude = 2 * np.abs(coefficients[k]) / n  # Amplitude berechnen
        phase = np.angle(coefficients[k])  # Phase berechnen
        result += amplitude * smp.cos(2 * 3.141592653589793 * frequencies[k] * r + phase) # Nur Realteil interessant also bleibt nur der cos(.)
        
    return result


# Fuktionswert an der Stelle r über eine Fourierreihe aus diskreten Inputwerten bestimmen
def FourierSeries(r, inputFunction):
    fftInputFunction = np.fft.fft(inputFunction) # Koeffizienten der fft berechnet mit der Funktion fft aus dem fft package
    return inverseFourier(r, fftInputFunction, np.fft.fftfreq(len(inputFunction)))

def dFourierSeries(r, inputFunction):
    fftInputFunction = np.fft.fft(inputFunction)
    coefficients = fftInputFunction
    frequencies = np.fft.fftfreq(len(inputFunction))
    
    n = len(coefficients)
    result = np.real(coefficients[0]) / n # konstanter Anteil
    
    for k in range(1, n//2):  # aus symmetriegründen nur die erste Hälfte der Freuqenzen zu betrachten
        amplitude = 2 * np.abs(coefficients[k]) / n  # Amplitude berechnen
        phase = np.angle(coefficients[k])  # Phase berechnen
        result += - amplitude * 2 * 3.141592653589793 * coefficients[k] * smp.sin(2 * 3.141592653589793 * frequencies[k] * r + phase) # wegen i als Vorfaktor durch die Ableitung ist für den Realteil nur -sin(.) interessant
        
    return result

# Funktion die das Anfangswertproblem definiert zum lösen mit odeitn
def dSdr(r, S):
    s_r, s_phi = S
    return [    1/r * s_phi   - calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   - 1/r * s_r,
                - 1/r * s_phi  + calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   + 1/r * s_r   + fourierFunctionNy_f(r) * ds_z   + dfourierFunctionNy_f(r) * (s_r + s_z)   + dfourierFunctionE_f(r) * 1/fourierFunctionE_f(r) * (-fourierFunctionNy_f(r) * s_r  + s_phi  + fourierFunctionNy_f(r) * s_z)   - (1 + fourierFunctionNy_f(r)) * radialForce_f(r)]


# input Function for solbvp with E(r) and ny(r) as Fourrier Series
def funcFourier(r, y):
    dSigmadr = 1/r * y[1]   - calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   - 1/r * y[0],
    dSigmaPfi = -1/r * y[1]   + calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   + 1/r * y[0]   + fourierFunctionNy_f(r) * ds_z   + dfourierFunctionNy_f(r) * (y[0] + s_z)   + dfourierFunctionE_f(r) * 1/fourierFunctionE_f(r) * (-fourierFunctionNy_f(r) * y[0]  + y[1]  + fourierFunctionNy_f(r) * s_z)   - (1 + fourierFunctionNy_f(r)) * radialForce_f(r)
    
    return np.vstack((dSigmadr,dSigmaPfi))


# definition of the boundary conditions for solving the bvp eacht element in the return array will be equal to zero for the solution
# with y[0] as s_r and y[1] as s_phi ans ya on the inner side and yb on the outer side
def bc(ya, yb):
    # return np.array([ya[0],ya[1]-4.2E08])
    return np.array([ya[0], yb[0]])
 
# input Function for solbvp with dE(r) and dny(r) calculated by the gradient
def funcGradient(r, y):
    #print("r: ",r[:10])
    #print("r[len(r)- 5:]", r[len(r)- 5:])
    #print(np.shape(r), "r2")
    # print(np.shape(y), "y")
    # #r = np.linspace(r_i, r_a, len(y[1]))
    # print(np.shape(r), "rlinspace")
    # print(np.shape(s_z))
    # print(s_z)
    ##################################################################################################################
    #lengthR = len(r)
    #E = np.ones(lengthR) * 100 * 10**9 * mSkalierung**(-1)  # E Modul in N/m^2 * mSaklierung**(-2)
    #E[:int(0.2* lengthR)] = 150 * 10**9 * mSkalierung**(-1)
    #E[len(E) - int(0.2* lengthR):] = 150 * 10**9 * mSkalierung**(-1)
    E = calcTanhValue(r, coefficientsE, materialEs)
    ny = calcTanhValue(r, coefficientsNy, materialNys)
    #ny = np.ones(lengthR) * 0.3 # Querkontraktionszahl
    #ny[len(E) - int(0.2* lengthR):] = 0.25
    #ny[:int(0.2* lengthR)] = 0.25
    
    dGradientE = np.gradient(E)
    dGradientNy = np.gradient(ny)

    dSigmaR = 1/r * y[1]   - calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   - 1/r * y[0],
    dSigmaPfi = (-1/r * y[1]   + calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   + 1/r * y[0]   
                   + ny * ds_z   + dGradientNy * (y[0] + s_z)   
                   + dGradientE * 1/E * (-ny * y[0]  + y[1]  + ny * s_z)   
                   - (1 + ny) * calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j )
    
    #print(np.shape(E[numberOfValues - len(r):]), "E")
    #print(np.shape(-1/r * y[1]   + calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   + 1/r * y[0]))
    #print(np.shape(+ ny * ds_z   + dGradientNy * (y[0] + s_z) ))
    #print(np.shape(+ dGradientE * 1/E * (-ny * y[0]  + y[1]  + ny * s_z)  ))
    #print(np.shape(- (1 + ny) * calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j ))
    # dSigmadr = 1/r * y[1]   - calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   - 1/r * y[0],
    # dSigmaPfidr = (-1/r * y[1]   + calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   + 1/r * y[0]   
    #                + ny[numberOfValues - len(r):] * ds_z   + dGradientNy[numberOfValues - len(r):] * (y[0] + s_z)   
    #                + dGradientE[numberOfValues - len(r):] * 1/E[numberOfValues - len(r):] * (-ny[numberOfValues - len(r):] * y[0]  + y[1]  + ny[numberOfValues - len(r):] * s_z)   
    #                - (1 + ny[numberOfValues - len(r):]) * calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j )
    
    return np.vstack((dSigmaR,dSigmaPfi))
# def funcGradient(r, y):
#     print("r: ",r[:10])
#     print("r[len(r)- 5:]", r[len(r)- 5:])
#     print(np.shape(r), "r2")
#     #print(np.shape(E[numberOfValues - len(r):]), "E")
#     #print(np.shape(-1/r * y[1]   + calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   + 1/r * y[0]))
#     #print(np.shape(+ ny * ds_z   + dGradientNy * (y[0] + s_z) ))
#     #print(np.shape(+ dGradientE * 1/E * (-ny * y[0]  + y[1]  + ny * s_z)  ))
#     #print(np.shape(- (1 + ny) * calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j ))
#     dSigmadr = 1/r * y[1]   - calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   - 1/r * y[0],
#     dSigmaPfidr = (-1/r * y[1]   + calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   + 1/r * y[0]   
#                    + ny[numberOfValues - len(r):] * ds_z   + dGradientNy[numberOfValues - len(r):] * (y[0] + s_z)   
#                    + dGradientE[numberOfValues - len(r):] * 1/E[numberOfValues - len(r):] * (-ny[numberOfValues - len(r):] * y[0]  + y[1]  + ny[numberOfValues - len(r):] * s_z)   
#                    - (1 + ny[numberOfValues - len(r):]) * calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j )
    
#     return np.vstack((dSigmadr,dSigmaPfidr))

# # Vektor containing s_r and ds_r/dr = h
# def dSdr(r, S):
#     s_r, h = S
#     return [    h,
#                 -1/r * 3 * h    - dradialForce_f(r)    - 1/r * (+fourierFunctionNy(r) * ds_z + dfourierFunctionNy_f(r) * (s_r + s_z) + dfourierFunctionE_f(r) * 1/fourierFunctionE(r) * (-fourierFunctionNy(r) * s_r + r * h + r*radialForce_f(r) + s_r - fourierFunctionNy(r) * s_z) - (2 + fourierFunctionNy_f(r)) * radialForce_f(r))] 
#


def calcDisplacement(r, s_r, s_phi, E, ny):
    #print(type(E), type(ny), type(s_r), type(s_phi), type(r))
    #print(np.size(E), np.size(ny), np.size(s_r), np.size(s_phi), np.size(r))
    u_r = ( 1/E * (-ny * s_r + s_phi) ) * r
    #print("ccccccccc")
    e_r = 1/E * (s_r - ny * s_phi)
    #e_r = r
    #print("test2")
    return([u_r, e_r])

####################################################################################################################
######  Hauptprogramm

## Berechnen der differenzierbaren Funktionen E(r) und Ny(r) und ihrer ersten Ableitungen nach r
#print("Berechne Fourrierreihen für E, ny und R und die entsprechenden Ableitungen")

# E(r)
x = smp.Symbol("x", real=True)
# fourierFunctionE = FourierSeries(x, E) # A function dependant on the "Symbol" x
# #print(FourierFunctionE)
# fourierFunctionE_f = smp.lambdify(x, fourierFunctionE) # Convert FourierFunctionEr to a numerical function for plotting
# #FourierFunctionEr = FourierSeries(np.linspace(0, len(r), len(r)*5), Er)


# # ny(r)
# fourierFunctionNy = FourierSeries(x, ny) # A function dependant on the "Symbol" x
# fourierFunctionNy_f = smp.lambdify(x, fourierFunctionNy) # Convert FourierFunctionEr to a numerical function for plotting
# #FourierFunctionNy = FourierSeries(np.linspace(0, len(r), len(r)*5), ny)


# # dE(r)/dr
# dfourierFunctionE = smp.diff(fourierFunctionE, x)
# dfourierFunctionE_f = smp.lambdify(x, dfourierFunctionE)
# dFormulaFourierFunctionE = dFourierSeries(x, E)
# dFormulaFourierFunctionE_f = smp.lambdify(x, dFormulaFourierFunctionE)

# # dny(r)/dr
# dfourierFunctionNy = smp.diff(fourierFunctionNy, x)
# dfourierFunctionNy_f = smp.lambdify(x, dfourierFunctionNy)
# dFormulaFourierFunctionNy = dFourierSeries(x, ny)
# dFormulaFourierFunctionNy_f = smp.lambdify(x, dFormulaFourierFunctionNy)

### Bestimmen der Spannungen in radiale Richtungen s_r(r) und in Umfangsrichtung s_phi(r)

# calculating the radial force, the numerical pendant, the first derivative and the numerical pendant
radialForce = calcBFeld(x, r_a, r_i, b_za, b_zi, b_0) * j
radialForce_f = smp.lambdify(x, radialForce)
dradialForce = smp.diff(radialForce, x)
dradialForce_f = smp.lambdify(x, dradialForce)

####################################################################################################################
####   lösen des DGL Systems 1. Ordnung bestehend aus (b*) für d(sigma_r)/dr und (a*) für d(sigma_phi)/dr
print("Lösen des DGL Systems 1. Ordnung:")

# ###   solving initial value problem with odeint
# print("mit odeint")
# # defining initial conditions for s_r(r_i) and s_phi(r_i)   (initial conditions are the values for r[0])
# S_0 = (s_r0, s_phi0)
# solOdeint = odeint(dSdr, y0=S_0, t=r, tfirst=True)


# s_rsolOdeint = solOdeint.T[0]
# s_phisolOdeint = solOdeint.T[1]


###   solving bvp with solve_bvp 
print("mit solvebvp")
##   mit Fourrierreiehn
print("für Fourrierreihen")

# Sguess = np.zeros((len(r), len(r)), dtype=float)
# yInitialGuess = np.zeros((2, r.size))

# solBvpFourier = solve_bvp(funcFourier, bc, r, yInitialGuess)

# s_rSolBvpFourier = solBvpFourier.sol(r)[0]
# s_phiSolBvpFourier = solBvpFourier.sol(r)[1]


## mit gradient
print("für die Recheckfunktionen")


yInitialGuess = np.zeros((2, r.size))

dGradientE = np.gradient(E)
dGradientNy = np.gradient(ny)
#print(np.shape(dGradientE))
#print(np.shape(E))
#print(np.shape(dGradientNy))
#print(np.shape(ny))
#print(np.shape(r), "r")
solBvpGradient = solve_bvp(funcGradient, bc, r, yInitialGuess, max_nodes=numberOfValues * 2)
print("solbvp results:\n solStatus (0 wenn alles geklappt hat): ", solBvpGradient.status, "\nsolmessage: ", solBvpGradient.message)

s_rSolBvpGradient = solBvpGradient.sol(r)[0]
s_phiSolBvpGradient = solBvpGradient.sol(r)[1]

# print("für die tanh Funktion")

# dGradientE = np.gradient(E)
# dGradientNy = np.gradient(ny)
# #print(np.shape(dGradientE))
# #print(np.shape(E))
# #print(np.shape(dGradientNy))
# #print(np.shape(ny))
# #print(np.shape(r), "r")
# solBvpGradient = solve_bvp(funcGradient, bc, r, yInitialGuess, max_nodes=10000)
# print("solbvp results:\n solStatus (0 wenn alles geklappt hat): ", solBvpGradient.status, "\nsolmessage: ", solBvpGradient.message)

# s_rSolBvpGradient = solBvpGradient.sol(r)[0]  * mSkalierung**(-1)  # für die Spannungen gilt kg/(s^2mm) = kg/(s^2 mm^2/mm) = (kg mm/s^2) / mm^2 = (kg m * 10-3 /s^2) / mm^2 = N/mm^2 * 10-3
# s_phiSolBvpGradient = solBvpGradient.sol(r)[1]  * mSkalierung**(-1)   # für die Spannungen gilt kg/(s^2mm) = kg/(s^2 mm^2/mm) = (kg mm/s^2) / mm^2 = (kg m * 10-3 /s^2) / mm^2 = N/mm^2 * 10-3


# ###################### lösen der Gleichung (e) also einer DGL 2. Ordnung durch überführen in ein DGL System mit DGLs 1. Ordnung


####################################################################################################################
### Berechnen der Verschiebungen
print("Berechnen der Verschiebungen")

u_rGradient = calcDisplacement(r, s_rSolBvpGradient, s_phiSolBvpGradient, E, ny)[0]
# u_rCaldwell = (calcDisplacement(r * mSkalierung*(-1),  calcStresses(r=r * mSkalierung**(-1), s_z0=0, s_ri=0, s_ra=0, nu=0.3, b_za=b_za, b_zi=b_zi, b_0=b_0, j=j * mSkalierung**(2))[0],
#                                                        calcStresses(r=r * mSkalierung**(-1), s_z0=0, s_ri=0, s_ra=0, nu=0.3, b_za=b_za, b_zi=b_zi, b_0=b_0, j=j * mSkalierung**(2))[1], 
#                                100 * 10**9, 0.3)[0] #########################################################################################################
#                 * mSkalierung**1)      # Berechnung von calc Displacement in SI Basiseinheiten und dann Umrechnung von m in mm
u_rCaldwell = 1000**(-1) * calcDisplacement(r * mSkalierung**(-1), mSkalierung**(-1) * calcStresses(r=r * mSkalierung**(-1), r_a=r_a * mSkalierung**(-1), r_i=r_i * mSkalierung**(-1), s_z0=s_z0, s_ri=s_ri, s_ra=s_ra, nu=ny[0], b_za=b_za, b_zi=b_zi, b_0=b_0, j=j * mSkalierung**2)[0],
                                                                   mSkalierung**(-1) * calcStresses(r=r * mSkalierung**(-1), r_a=r_a * mSkalierung**(-1), r_i=r_i * mSkalierung**(-1), s_z0=s_z0, s_ri=s_ri, s_ra=s_ra, nu=ny[0], b_za=b_za, b_zi=b_zi, b_0=b_0, j=j * mSkalierung**2)[1],
                                E[0], ny[0])[0]

##### Plots
#rTest = np.linspace(r_i, r_a, 5000)
# sArray = calcStresses(rTest, 0, 0, 0, 0.3, b_za, b_zi, b_0, j)

# ax1 = plt.subplot(2,1,1)
# ax1.plot(rTest, sArray[0])
# ax2 = plt.subplot(2,1,2)
# ax2.plot(rTest, sArray[1])
# plt.title("wirklich jonas")
# #ax1.xlabel("testtesttest")
# plt.show()



# E(r)
plt.subplot(anzahlRowPlots, anzahlColumnPlots, 1)
plt.plot(r , E, label="vorgegebenes E(r)")
plt.xlabel(f"Radius in m E{int(np.log10(mSkalierung))}")
plt.ylabel(f"E-Modul in N/(m E{int(np.log10(mSkalierung))})^2")

#plt.plot(np.linspace(r_i, r_a, numberOfValues*5), 
#         fourierFunctionE_f(x=np.linspace(0, len(r), len(r) * 5)), label='E(r) mit Fourrierreihe genähert', linestyle='--')
#plt.legend()


# ny(r)
plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 1)
plt.plot(r , ny, label="vorgegebenes ny(r)")
plt.xlabel(f"Radius in m E{int(np.log10(mSkalierung))}")
plt.ylabel("ny in 1")

# plt.plot(np.linspace(r_i, r_a, numberOfValues*5), 
#          fourierFunctionNy_f(x=np.linspace(0, len(r), len(r) * 5)), label='ny(r) mit Fourrierreihe genähert', linestyle='--')
# plt.legend()


plt.subplot(anzahlRowPlots, anzahlColumnPlots, 2)
# plt.plot(r, s_rsolOdeint, "--", label="s_r berechnet als IVP")
plt.xlabel(f"Radius in m E{int(np.log10(mSkalierung))}")
plt.ylabel(f"s_r in N/(m E{int(np.log10(mSkalierung))})^2")
plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 2)
# plt.plot(r, s_phisolOdeint, "--", label="s_phi berechnet als IVP")
plt.xlabel(f"Radius in m E{int(np.log10(mSkalierung))}")
plt.ylabel(f"s_phi in N/(m E{int(np.log10(mSkalierung))})^2")

# print(mSkalierung**(-2) * calcStresses(r=r * mSkalierung**(-1), s_z0=s_z0, s_ri=s_ri, s_ra=s_ra, nu=ny[0], b_za=b_za, b_zi=b_zi, b_0=b_0, j=j * mSkalierung**2)[0])
#print(mSkalierung**(-2) * calcStresses(r=r * mSkalierung**(-1), s_z0=s_z0, s_ri=s_ri, s_ra=s_ra, nu=ny[0], b_za=b_za, b_zi=b_zi, b_0=b_0, j=j * mSkalierung**2)[0])


plt.subplot(anzahlRowPlots, anzahlColumnPlots, 2)
#plt.plot(r, s_rSolBvpFourier, label="s_r berechnet als BVP mit Fourierreihe")
plt.plot(r, mSkalierung**(-2) * calcStresses(r=r * mSkalierung**(-1), r_a=r_a * mSkalierung**(-1), r_i=r_i * mSkalierung**(-1), s_z0=s_z0, s_ra=s_ra, s_ri=s_ri, nu=ny[0], b_za=b_za, b_zi=b_zi, b_0=b_0, j=j * mSkalierung**2) [0], label="s_r nach Caldwell")
# plt.plot(rTest, calcStresses(rTest,r_a,r_i, 0, 0, 0, 0.3, b_za, b_zi, b_0, j)[0], label="s_r nach Caldwell")

plt.plot(r, s_rSolBvpGradient * mSkalierung**(-1), "--", label="s_r berechnet als BVP mit Differenzenquotienten") 
 # für die Spannungen gilt kg/(s^2mm) = kg/(s^2 mm^2/mm) = (kg mm/s^2) / mm^2 = (kg m * 10-3 /s^2) / mm^2 = N/mm^2 * 10-3
plt.legend()

plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 2)
#plt.plot(r, s_phiSolBvpFourier, label="s_phi berechnet als BVP mit Fourierreihe")
plt.plot(r, mSkalierung**(-2) * calcStresses(r=r * mSkalierung**(-1), r_a=r_a * mSkalierung**(-1), r_i=r_i * mSkalierung**(-1), s_z0=s_z0, s_ri=s_ri, s_ra=s_ra, nu=ny[0], b_za=b_za, b_zi=b_zi, b_0=b_0, j=j * mSkalierung**2)[1], label="s_phi nach Caldwell")
plt.plot(r, s_phiSolBvpGradient * mSkalierung**(-1), "--", label="s_phi berechnet als BVP mit Differenzenquotienten") 
# für die Spannungen gilt kg/(s^2mm) = kg/(s^2 mm^2/mm) = (kg mm/s^2) / mm^2 = (kg m * 10-3 /s^2) / mm^2 = N/mm^2 * 10-3
plt.legend()


plt.subplot(anzahlRowPlots, anzahlColumnPlots, 3)
#plt.plot(r, dfourierFunctionE_f(r), label="dFourierE")
#plt.plot(r, dfourierFunctionNy_f(r) * 10**6, label="dFourierNy * 10^6")
#plt.plot(r, dFormulaFourierFunctionE_f(r), label="dFormulaFourierE")
plt.xlabel(f"Radius in m E{int(np.log10(mSkalierung))}")
plt.ylabel(f"Änderung für dE in N/mm^3 für dny in 1/m E{int(np.log10(mSkalierung))} E2")
#plt.plot(r, radialForce_f(r), label="R")
plt.legend()
#plt.plot(r, dradialForce_f(r))


plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 3)
plt.plot(r, dGradientE * 10**3, label="dGradientE * 10^3")
plt.plot(r, dGradientNy * 10**6, label="dGradientNy * 10^6")
plt.plot(r, calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j * 10, label="R * 10")
plt.xlabel(f"Radius in m E{int(np.log10(mSkalierung))}")
plt.ylabel(f"Änderung für dE in N/(m E{int(np.log10(mSkalierung))})^3 für dny in 1/m E{int(np.log10(mSkalierung))}")
plt.legend()


plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 4)
#u_rFourierFunction = calcDisplacement(r, s_rSolBvpFourier, s_phiSolBvpFourier, fourierFunctionE_f(r), fourierFunctionNy_f(r))[0]
#plt.plot(r, u_rFourierFunction, color="orange", label="u_r berechnet über BVP mit Fourierreihe",)
plt.plot(r, u_rGradient, label="u_r berechnet über BVP mit Differenzenquotient")
#plt.plot(r, u_rCaldwell, "b", label="u_r berechnet über Caldwell")
plt.plot(r, u_rCaldwell, "b", label="u_r berechnet über Caldwell")
#plt.plot(r, u_rCaldwell)
plt.xlabel(f"Radius in m E{int(np.log10(mSkalierung))}")
plt.ylabel(f"u_r in m E{int(np.log10(mSkalierung))}")
plt.legend()


# Save plot
fig = plt.gcf()
fig.set_size_inches(20, 12)
currentTime = datetime.now()
pictureName = f"Bilder\Graphen{currentTime.year}-{currentTime.month:02d}-{currentTime.day:02d}-{currentTime.hour:02d}-{currentTime.minute:02d}"
plt.savefig(pictureName, dpi=500)

print("Fertig")
plt.show()