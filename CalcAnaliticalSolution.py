import matplotlib.pyplot as plt
import numpy as np
import sympy as smp
import math as m
from scipy.sparse import *
import scipy.sparse as sparse
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import spsolve
from scipy.linalg import det
from scipy.integrate import odeint
from scipy.integrate import solve_bvp
from functions.Hoopstress_Calc02Copy import calcStresses
from functions.calcMaterialFunction import calcMaterialTanhCoefficients
from functions.calcMaterialFunction import calcTanhValue
from functions.calcMaterialFunction import calcTanhDerivativeValue
from datetime import datetime # Nur für die Bennenung der Grafiken

#####   Definierte Werte
mSkalierung = 1        # Skalierungsfaktor für die Einheit Meter gleich 1, für mm gleich 10*(3), etc. 
windingDivisor = 60        # Es wird für 600/windingDivisor Windungen gerechnet
numberOfValues = int(300000 / windingDivisor) # Anzahl der Werte des r Vektors
solbvpMaxNodes = numberOfValues * 3
totalWindings = int(600 /windingDivisor)

#   Subplots
anzahlRowPlots = 2
anzahlColumnPlots = 4

# #   Spannungsrandbedingungen

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
materialEs = [E_con, E_cow, E_ins]

ny_con = 0.35 #* 10**9 * mSkalierung**(-1) # Possion's Ratio Conductor
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


# # Funktion zur Berechnung des B-Felds an der Stelle r
# def calcBFeld(r, r_a, r_i, b_za, b_zi, b_0):
#     return ((b_za - b_zi)/(r_a - r_i))  *  r + b_0

# coefficientsE = calcMaterialTanhCoefficients(r_i , r_a, t, materialWidths, materialEs, slope=1000, scale=720)
# E = calcTanhValue(r, coefficientsE, materialEs)
# print(E)
# #E = np.ones(numberOfValues) * 100 * 10**9 * mSkalierung**(-1) # E Modul in  N/m^2 E-3 = kg m/(s^2 m^2) E-3 = kg/(s^2 m) E-3 = kg/(s^2 mm)
# #E[:int(0.2* numberOfValues)] = 150 * 10**9 * mSkalierung**(-1)
# #E[len(E) - int(0.2* numberOfValues):] = 150 * 10**9 * mSkalierung**(-1)


# coefficientsNy = calcMaterialTanhCoefficients(r_i, r_a, t, materialWidths, materialNys, slope=1000, scale=720)
# ny = calcTanhValue(r, coefficientsNy, materialNys)
# #ny = np.ones(numberOfValues) * 0.3 # Querkontraktionszahl
# #ny[len(E) - int(0.2* numberOfValues):] = 0.25
# #ny[:int(0.2* numberOfValues)] = 0.25


def calcIntegral(n, r, r_Start, r_Exterior, r_Center, j, b_za, b_zi, b_0):
    return j * (  (b_za - b_zi)/(r_Exterior - r_Center) * (1 / (n+2) ) * (r**(n+2) - r_Start**(n+2)) 
                 +  b_0                     * (1/  (n+1) ) * (r**(n+1) - r_Start**(n+1)))
    

def calcAnaliticalSolution(rDomains, rCenter, rExterior, s_rCenter, s_rExterior, s_zBegin, windings, nyArray, j, b_za, b_zi, b_0):
    # rDomains = [ [r_1i] , [r_2i] , ...., [r_ni]] with np.arrays; endpoints of domain not in array
    ### Berechne s_z
    rArray = rDomains.flatten()
    s_z = s_zBegin  # from eq (4) ################################################################
    print(type(rDomains))
    rEnds = [item[0] for item in rDomains]
    rEnds.pop(0)
    rEnds.append(rExterior)
    #print(rEnds)
    # Equation (14) with sigma_z=const.
    cB = lambda r, r1 : 1/2 * (1 - r1**2/r**2)
    cA = lambda r, r1 : 1/2 * (1 + r1**2/r**2)
    c = lambda r, r1, ny : - 1/2 * (1+ny) * calcIntegral(n=0, r=r, r_Start=r1, r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0) - 1/r**2 * (1 - (1+ny)/2 ) * calcIntegral(n=2, r=r, r_Start=r1, r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
    row = []
    col = []  # column
    data = []
    constant = []
    nyTestDummy = []
    for i, rJump in enumerate(rEnds):
        if i == 0:
            riStart = rCenter
            constant.append(cA(rJump, riStart) * s_rCenter + c(rJump, riStart, nyRespectR( (riStart+rJump)/2 , rArray, nyArray)) ) 
        else:
            riStart = rEnds[i - 1]
            constant.append(c(rJump, riStart, nyRespectR( (riStart+rJump)/2 , rArray, nyArray)))
        row.extend([i, i, i])
        icol = 3*i - 1
        col.extend([icol, icol +1, icol +3])
        data.extend([-cA(rJump, riStart), -cB(rJump, riStart), 1])
        #print(nyRespectR((riStart+rJump)/2 , rArray, nyArray))
        nyTestDummy.append(nyRespectR((riStart+rJump)/2 , rArray, nyArray))
    print("first part of nyTestDummy ", nyTestDummy[:20])
    print("last part of nyTestDummy ", nyTestDummy[-20:])
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
    length = len(data)
    nyTestDummy2 = []
    for i, rJump in enumerate(rEnds):
        if i == 0:
            riStart = rCenter                                                                               
            constant.append(s_rCenter - (1+nyRespectR( (riStart+rJump)/2 , rArray, nyArray)) * calcIntegral(n=0, r=rJump, r_Start=riStart, r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)) 
        else:
            riStart = rEnds[i - 1]
            constant.append(- (1+nyRespectR( (riStart+rJump)/2 , rArray, nyArray)) * calcIntegral(n=0, r=rJump, r_Start=riStart, r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0))
        row.extend([windings*len(materialWidths) + i,  windings*len(materialWidths) + i,  windings*len(materialWidths) + i,  windings*len(materialWidths) + i])
        icol = 3*i - 1
        col.extend([icol, icol +1, icol +2, icol +3])
        data.extend([-1, -1, 1, 1])
        #print(nyRespectR((riStart+rJump)/2 , rArray, nyArray))
        nyTestDummy2.append(nyRespectR((riStart+rJump)/2 , rArray, nyArray))
        if nyTestDummy[i] != nyRespectR((riStart+rJump)/2 , rArray, nyArray):
            print("Error !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    # remove first and last element of the newly added values because they are know (sigma_r(r_center) and sigma_r(r_exterior))
    print("first part of nyTestDummy2", nyTestDummy2[:20])
    print("last part of nyTestDummy2", nyTestDummy2[-20:])
    row.pop(length)    # not (length - 1) because the first element added after length was assigned shall be poped ##########################################
    col.pop(length)
    data.pop(length)
    row.pop()
    col.pop()
    data.pop()
    # add -sigma_r(r_exterior) to the constant vektor
    constant[-1] += -s_rExterior
    
    # u(r_i-1End) = u_(r_iStart) with sigma_z=const.
    length = len(data)
    rEndsWithoutLastEnd = rEnds[:] # [:] needed to copy the list and not assigning a reference to the original list
    rEndsWithoutLastEnd.pop()
    print("len(rEnds) len(rEndsWithoutLastEnd)", len(rEnds), len(rEndsWithoutLastEnd))
    for i, rJump in enumerate(rEndsWithoutLastEnd):       # iterating over the material transitions ##########################
        if i == 0:
            riStart = rCenter
 #       elif:
        else:
            riStart = rEnds[i - 1]
        row.extend([2 * windings*len(materialWidths) + i, 2 * windings*len(materialWidths) + i, 2 * windings*len(materialWidths) + i])
        icol = 3*i + 1
        col.extend([icol, icol +1, icol +2])
        data.extend([1, -nyRespectR( (rJump+riStart)/2 , rArray, nyArray) + nyRespectR( (rEnds[i+1]+rJump)/2 , rArray, nyArray), -1])###################################################
        #print(nyRespectR((riStart+rJump)/2 , rArray, nyArray))
        if nyTestDummy[i] != nyRespectR((riStart+rJump)/2 , rArray, nyArray): #nur zum prüfen von ny
            print("Error !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    
    constant.extend([0] * (len(rEnds) - 1))
    
    # for i, rJump in enumerate(rEnds):
    #     if i != 0:                      # iterating over the material transitions ########################## Better for i, rJump in enumerate(rEnds.pop(0)): i += 1
    #         riStart = rEnds[i - 1]
    #         i = i - 1
    #         row.extend([windings*len(materialWidths) - 1 + i, windings*len(materialWidths) - 1 + i, windings*len(materialWidths) - 1 + i])
    #         icol = 3*i + 1
    #         col.extend([icol, icol +1, icol +2])
    #         data.extend([1, -nyRespectR( (rJump-rEnds[i])/2 ) + nyRespectR( (rEnds[i+1]-rJump)/2 ), -1])###################################################
    # constant.extend([0] * (len(rEnds) - 1))
    
    # calculate the sparse matrix
    A = csr_matrix((data, (row, col)), shape=(3*windings*len(materialWidths) -1, 3*windings*len(materialWidths) -1), dtype=float)
    #print(A.toarray())
    print(A.shape)
    print(type(A))
    #print(constant)
    print(len(constant))
    print(type(constant))
    
    # calculate all Streses at the material transitions
    #print("A in sparse: ", A[windings*len(materialWidths) - 5:])
    #print("A in dense: ", A.todense())
    print(type(A.todense()))
    print("Determinant: ", det(A.todense()))
    
    
    nonzero_indices = A.nonzero()
    nonzero_values = A.data
    #print(len(data), len(nonzero_values))
    
    plt.scatter(col, row, c=nonzero_values, cmap="YlGnBu",)
    plt.grid(True)
    plt.colorbar()
    plt.ylabel("Zeile")
    plt.gca().invert_yaxis()
    plt.xlabel("Spalte")
    plt.show()
    
    
    Ainv = inv(A)
    #Ainv = spsolve(A, identity) # identity is from scipy.sparse
    result = Ainv.dot(constant)
    print(result)
        
    ###Berechnen der Spannungen im innersten Torus
    #r = np.linspace(rCenter, rEnds[0], rValuesPerMaterial)
    s_rArray = (  (1/2) * result[0] * (1 - rCenter**2/rDomains[0]**2)                               
                - ((1 + nyRespectR(rCenter, rArray, nyArray)) / 2) * calcIntegral(n=0, r=rDomains[0], r_Start=rDomains[0][0], r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                -  (1/rDomains[0]**2) * (1 - (1 + nyRespectR(rCenter, rArray, nyArray)) /2 ) * calcIntegral(n=2, r=rDomains[0], r_Start=rDomains[0][0], r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                +  (s_rCenter  *  (1/2)  *  (1 + (rCenter**2/rDomains[0]**2))))
    #print("s_rArray: ", s_rArray[0:5])
    
    
    ### Berechne s_phi
    s_phiArray = ( (result[0] - nyRespectR(rCenter, rArray, nyArray) * s_zBegin)  
                 +  nyRespectR(rCenter, rArray, nyArray) * s_z  
                 -  (1 + nyRespectR(rCenter, rArray, nyArray)) * calcIntegral(n=0, r=rDomains[0], r_Start=rDomains[0][0], r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                 + s_rCenter
                 - s_rArray )
    
    for winding in range(len(rDomains) - 1):
        #rNew = np.linspace()######################################################
        winding += 1
        # print(result[3*winding], rEnds[winding - 1], len(rDomains[winding]))
        # print(len(- ((1 + nyRespectR( (rEnds[winding]-rEnds[winding-1]/2), rArray, nyArray)) / 2) * calcIntegral(0, rDomains[winding], rDomains[winding][0], rExterior, rCenter, j, b_za, b_zi, b_0)))
        # print(len(-  (1/rDomains[winding]**2) * (1 - (1 + nyRespectR( (rEnds[winding]-rEnds[winding-1]/2), rArray, nyArray))/2 ) * calcIntegral(2, rDomains[winding], rDomains[winding][0], rExterior, rCenter, j, b_za, b_zi, b_0)))
        # print(len(+  (result[3*winding - 1]  *  (1/2)  *  (1 + (rEnds[winding - 1]**2/r**2)))))
        # print(+  (result[3*winding - 1]  *  (1/2)  ))
        # print((1 + (rEnds[winding - 1])))
        # print(len(r**2))
        
        s_rArrayNew = (  (1/2) * result[3*winding] * (1 - rEnds[winding - 1]**2/rDomains[winding]**2) 
                    - ((1 + nyRespectR( (rEnds[winding]+rEnds[winding-1])/2, rArray, nyArray)) / 2) * calcIntegral(n=0, r=rDomains[winding], r_Start=rDomains[winding][0], r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                    -  (1/rDomains[winding]**2) * (1 - (1 + nyRespectR( (rEnds[winding]+rEnds[winding-1])/2, rArray, nyArray)) /2 ) * calcIntegral(n=2, r=rDomains[winding], r_Start=rDomains[winding][0], r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                    +  (result[3*winding - 1]  *  (1/2)  *  (1 + (rEnds[winding - 1]**2/rDomains[winding]**2))))
        
        s_phiArrayNew = ( (result[3*winding] - nyRespectR( (rEnds[winding]+rEnds[winding-1])/2 , rArray, nyArray) * s_zBegin)  
                        +  nyRespectR( (rEnds[winding]+rEnds[winding-1])/2, rArray, nyArray) * s_z  
                        -  (1 + nyRespectR( (rEnds[winding]+rEnds[winding-1])/2, rArray, nyArray)) * calcIntegral(n=0, r=rDomains[winding], r_Start=rDomains[winding][0], r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                        + result[3*winding - 1]
                        - s_rArrayNew )
        
        s_rArray = np.append(s_rArray, s_rArrayNew)
        s_phiArray = np.append(s_phiArray, s_phiArrayNew)
        #print(winding)
    
    
    #print(nyRespectR(r_i, rArray, nyArray))
    #print(r_i)
    return s_rArray, s_phiArray
    
def calcDomains(rCenter, materialWidths, numberOfWindings, numberOfValuesR):
    '''Calculates an np.array containing arrays with discrete values coresponding to the material domains; only the last domain contains the end value'''
    thisSpot = rCenter
    results = []
    for i in range(numberOfWindings):
        for width in materialWidths:
            results.append(np.linspace(thisSpot, thisSpot + width, int(numberOfValuesR/(numberOfWindings * len(materialWidths))) ))
            thisSpot += width
    return np.array(results)

def calcNyWithDomains(rDomains, materialNys):
    ny = np.ones_like(rDomains)
    for i, nyValue in enumerate(materialNys):
        ny[i::len(materialNys)] = ny[i::len(materialNys)] * nyValue
    print(len(ny))
    return ny.flatten()

def nyRespectR(r, rArray, nyArray):
    return np.interp(r, rArray, nyArray)


r = np.linspace(r_i, r_a, 500000)
rDom = calcDomains(r_i, materialWidths, totalWindings, len(r))
Ny = calcNyWithDomains(rDom, materialNys)
s_r, s_phi = calcAnaliticalSolution(rDomains=rDom, rCenter=r_i, rExterior=r_a, s_rCenter=0, s_rExterior=0, s_zBegin=0, windings=totalWindings, nyArray=Ny, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
plt.subplot(2,1,1)
plt.grid(True)
plt.plot(rDom.flatten(), s_r)
plt.subplot(2,1,2)
plt.grid(True)
plt.plot(rDom.flatten(), s_phi)
plt.show()