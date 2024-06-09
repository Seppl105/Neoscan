import matplotlib.pyplot as plt
import numpy as np
import sympy as smp
import math as m
from scipy.sparse import *
from scipy.sparse.linalg import inv
from scipy.integrate import odeint
from scipy.integrate import solve_bvp
from functions.Hoopstress_Calc02Copy import calcStresses
from functions.calcMaterialFunction import calcMaterialTanhCoefficients
from functions.calcMaterialFunction import calcTanhValue
from functions.calcMaterialFunction import calcTanhDerivativeValue
from datetime import datetime # Nur für die Bennenung der Grafiken

#####   Definierte Werte
mSkalierung = 10**(3)        # Skalierungsfaktor für die Einheit Meter gleich 1, für mm gleich 10*(3), etc. 
windingDivisor = 60        # Es wird für 600/windingDivisor Windungen gerechnet
numberOfValues = int(300000 / windingDivisor) # Anzahl der Werte des r Vektors
solbvpMaxNodes = numberOfValues * 3

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


coefficientsNy = calcMaterialTanhCoefficients(r_i, r_a, t, materialWidths, materialNys, slope=1000, scale=720)
ny = calcTanhValue(r, coefficientsNy, materialNys)
#ny = np.ones(numberOfValues) * 0.3 # Querkontraktionszahl
#ny[len(E) - int(0.2* numberOfValues):] = 0.25
#ny[:int(0.2* numberOfValues)] = 0.25


def calcIntegral(n, r, r_Start, r_Exterior, r_Center, j, b_za, b_zi, b_0):
    return j * (  (b_za - b_zi)/(r_Exterior - r_Center) * (1 / (n+2) ) * (r**(n+2) - r_Start**(n+2)) 
                 +  b_0                     * (1/  (n+1) ) * (r**(n+1) - r_Start**(n+1)))
    

def calcAnaliticalSolution(r, rEnds, rCenter, rExterior, s_rCenter, s_rExterior, windings, materialNys):
     # Equation (14) with sigma_z=const.
    cB = lambda r, r1 : 1/2 * (1 - r1**2/r**2)
    cA = lambda r, r1 : 1/2 * (1 + r1**2/r**2)
    c = lambda r, r1, ny : - 1/2 * (1+ny) * calcIntegral(0, r, r1, rExterior, rCenter, j, b_za, b_zi, b_0) - 1/r**2 * (1 - (1+ny) / 2) * calcIntegral(0, r, r1, rExterior, rCenter, j, b_za, b_zi, b_0)
    row = []
    col = []  # column
    data = []
    constant = []
    for i, rJump in enumerate(rEnds):
        if i == 0:
            riStart = 0
            constant.append(cA(rJump, riStart) * s_rCenter + c(rJump, riStart, ny( (riStart+rJump)/2 )) ) 
        else:
            riStart = rJump - rEnds[i - 1]
            constant.append(c(rJump, riStart, ny( (riStart+rJump)/2 )))
        row.extend([i, i, i])
        icol = 3*i - 1
        col.extend([icol, icol +1, icol +2])
        data.extend([-cA(rJump, riStart), -cB(rJump, riStart), 1])
    # remove first and last element because they are know (sigma_r(r_center) and sigma_r(r_exterior))
    row.pop(0)
    col.pop(0)
    data.pop(0)
    row.pop()
    col.pop()
    data.pop()
    # add -sigma_r(r_exterior) to the constant vektor
    constant[-1] += -s_rExterior
    
    # Equation (13) with sigma_z=const.
    len = len(data)
    for i, rJump in enumerate(rEnds):
        if i == 0:
            riStart = 0
            constant.append(-s_rCenter - (1+ny( (riStart+rJump)/2 )) * calcIntegral(0, r, riStart, rExterior, rCenter, j, b_za, b_zi, b_0)) 
        else:
            riStart = rJump - rEnds[i - 1]
            constant.append(- (1+ny( (riStart+rJump)/2 )) * calcIntegral(0, r, riStart, rExterior, rCenter, j, b_za, b_zi, b_0))
        row.extend([len + i, len + i, len + i, len + i])
        icol = 3*i - 1
        col.extend([icol, icol +1, icol +2, icol +3])
        data.extend([-1, -1, 1, 1])
    # remove first and last element of the newly added values because they are know (sigma_r(r_center) and sigma_r(r_exterior))
    row.pop(len)###########################################
    col.pop(len)
    data.pop(len)
    row.pop()
    col.pop()
    data.pop()
    # add -sigma_r(r_exterior) to the constant vektor
    constant[-1] += -s_rExterior
    
    # u(r_i-1End) = u_(r_iStart) with sigma_z=const.
    len = len(data)
    for i, rJump in enumerate(rEnds):
        if i != 0:                      # only iterating over the material transitions
            riStart = rJump - rEnds[i - 1]
            i = i - 1
            row.extend([len + i, len + i, len + i])
            icol = 3*i + 1
            col.extend([icol, icol +1, icol +2])
            data.extend([1, -ny( (rJump-rEnds[i])/2 ) + ny( (rEnds[i+1]-rJump)/2 ), -1])###################################################
    constant.extend([0] * (len(rEnds) - 1))
    
    # calculate the sparse matrix
    A = csr_matrix((data, (row, col)), shape=(3*windings -1, 3*windings -1), dtype=float)
    print(A.toarray())
    print(A.shape())
    print(constant)
    print(constant.shape())
    
    # calculate all Streses at the material transitions
    Ainv = inv(A)
    result = Ainv.dot(constant)
    
        
    ###Berechnen der Spannungen im innersten Torus
    r = np.linspace(rCenter, rEnds[0], valuesiiiiiiiiiiiiiiiiiiiiiiii)
    s_rArray = (  (1/2) * result[0] * (1 - rCenter**2/r**2) 
                - ((1 + ny[0]) / 2) * calcIntegral(0, r, rExterior, rCenter, j, b_za, b_zi, b_0)
                -  (1/r**2) * (1 - (1 + ny[0])/2 ) * calcIntegral(2, r, rExterior, rCenter, j, b_za, b_zi, b_0)
                +  (s_rCenter  *  (1/2)  *  (1 + (rCenter**2/r**2))))
    #print("s_rArray: ", s_rArray[0:5])
    
    
    ### Berechne s_phi
    s_phiArray = ( (result[0] - ny[0] * s_z0)  
                 +  ny[0] * s_z  
                 -  (ny[0] + 1) * calcIntegral(0, r, rExterior, rCenter, j, b_za, b_zi, b_0)
                 + s_rCenter
                 - s_rArray )
    
    for 
    
calcAnaliticalSolution(r=r, rEnds=r_, rCenter=r_i, rExterior=r_a, windings=10, materialNys=))