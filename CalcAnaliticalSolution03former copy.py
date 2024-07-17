import matplotlib.pyplot as plt
import numpy as np
#import sympy as smp
#import math as m
from scipy.sparse import *
import scipy.sparse as sparse
from scipy.sparse.linalg import inv
#from scipy.sparse.linalg import spsolve
from scipy.linalg import det
#from scipy.integrate import odeint
#from scipy.integrate import solve_bvp
from functions.Hoopstress_Calc02Copy import calcStresses
#from functions.calcMaterialFunction import calcMaterialTanhCoefficients
#from functions.calcMaterialFunction import calcTanhValue
#from functions.calcMaterialFunction import calcTanhDerivativeValue
from datetime import datetime # Nur für die Bennenung der Grafiken
#from Hoopstress_Calc06Maze import calcDisplacement
#from Hoopstress_Calc06Maze import plotDeviation



#####   Definierte Werte
mSkalierung = 1        # Skalierungsfaktor für die Einheit Meter gleich 1, für mm gleich 10*(3), etc. 
windingDivisor = 60       # Es wird für 600/windingDivisor Windungen gerechnet
numberOfValues = int(3000000 / windingDivisor) # Anzahl der Werte des r Vektors
#solbvpMaxNodes = numberOfValues * 3
totalWindings = int(600/windingDivisor)
#lenRPlot = 5000000

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

#t = 0.36 * 10**(-3) * mSkalierung     # [m] Dicke einer Wicklung (Conductor, Cowinding, Insulation)
t_con = 0.12 * 10**(-3) * mSkalierung # [m] Dicke des Bandleiters
t_cow = 0.23 * 10**(-3) * mSkalierung # [m] Dicke des Cowindings
t_ins = 0.01 * 10**(-3) * mSkalierung # [m] Dicke der Isolation
t = t_con + t_cow + t_ins             # [m] Dicke einer Wicklung (Conductor, Cowinding, Insulation)
materialWidths = [t_con, t_cow, t_ins]

E_con = 280 * 10**9 * mSkalierung**(-1) # E-Modul Conductor
E_cow = 300 * 10**9 * mSkalierung**(-1) # E-Modul Cowinding
E_ins = 200 * 10**9 * mSkalierung**(-1) # E-Modul Insulation
materialEs = [E_con, E_cow, E_ins]

ny_con = 0.35 #* 10**9 * mSkalierung**(-1) # Possion's Ratio Conductor
ny_cow = 0.3 # Possion's Ratio Cowinding
ny_ins = 0.4 # Possion's Ratio Insulation
materialNys = [ny_con, ny_cow, ny_ins]

#confinement
t_confinement = -1#0.0005       # if set to -1 the calculation will run without a confinement
E_confinement = E_con
ny_confinement = ny_con

#   Konstanten des Versuchsaufbau

r_i = 0.430 * mSkalierung               # [m] innerer radius
r_a = 0.646 * mSkalierung               # [m] äußerer radius ohne confinement
if r_a != (r_i + t * ( ((r_a - r_i) / t) / windingDivisor)):
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n r_a und t Werte mit der Anzahl der Windungen ergeben nicht dev vorgegegebenen äußeren Radius\n vorgegebener Wert für r_a: ", r_a, "\n mit Dicken berechneter Wert ((r_i + t * ( ((r_a - r_i) / t) / windingDivisor)): ", (r_i + t * ( ((r_a - r_i) / t) / windingDivisor)), "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
r_a = (r_i + t * ( ((r_a - r_i) / t) / windingDivisor))


#r = np.linspace(r_i,r_a, numberOfValues) # array mit diskreten Radien

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
    

def calcAnaliticalSolution(rDomains, rCenter, rExterior, s_rCenter, s_rOuter, s_zBegin, windings, nyArray, EArray, materialWidths, j, b_za, b_zi, b_0):
    
    ### Berechne s_z
    s_z = s_zBegin  # from eq (4) ################################################################
    
    # rDomains = [ [r_1i] , [r_2i] , ...., [r_ni]] with np.arrays; endpoints of domain not in array
    rArray = np.concatenate(rDomains, axis=None).ravel() ################# ravel() may be needless occures multiple times
    #rArray = rDomains.flatten()
    print(type(rDomains))
    rEnds = [item[0] for item in rDomains]
    rEnds.pop(0)
    rEnds.append(rArray[-1]) #rEnds.append(rExterior)
    
    # if confinement is defined some parts of the calculation have to be adadpted
    # wheater the calculation has to be adapted is controled with the variable "confinementDefined"
    jConfinement = 0
    confinementDefined = 0
    if len(rDomains) != windings*len(materialWidths):  #rExterior != rArray[-1]:
        confinementDefined = 1
        print("confinement ist definiert mit rArray[-1]: ", rArray[-1], "und rExterior: ", rExterior)
    #print(rEnds)
    
    # lists to assemble the sparse matrix
    row = []
    col = []  # column
    data = []
    constant = [] # list of constants of the system of equations
    nyTestDummy = []
    
    
    # first set of equations is assembled (one for each intervall of r) derived from eq. (14CA) with sigma_z=const.
    cB = lambda r, r1 : 1/2 * (1 - r1**2/r**2)
    cA = lambda r, r1 : 1/2 * (1 + r1**2/r**2)
    c = lambda r, r1, ny, jC: - 1/2 * (1+ny) * calcIntegral(n=0, r=r, r_Start=r1, r_Exterior=rExterior, r_Center=rCenter, j=jC, b_za=b_za, b_zi=b_zi, b_0=b_0) - 1/r**2 * (1 - (1+ny)/2 ) * calcIntegral(n=2, r=r, r_Start=r1, r_Exterior=rExterior, r_Center=rCenter, j=jC, b_za=b_za, b_zi=b_zi, b_0=b_0)

    for i, rJump in enumerate(rEnds):
        if i == 0:
            riStart = rCenter
            constant.append(cA(rJump, riStart) * s_rCenter + c(rJump, riStart, nyRespectR( (riStart+rJump)/2 , rArray, nyArray), j)) 
        else:
            riStart = rEnds[i - 1]
            if confinementDefined == 1 and rJump == rEnds[-1]:
                jNew = jConfinement
                print("NUR EINMAL STEHT DAS HIER DA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            else:
                jNew = j
            constant.append(c(rJump, riStart, nyRespectR( (riStart+rJump)/2 , rArray, nyArray), jNew))
        row.extend([i, i, i])
        icol = 3*i - 1
        col.extend([icol, icol +1, icol +3])
        data.extend([-cA(rJump, riStart), -cB(rJump, riStart), 1])
        #print(nyRespectR((riStart+rJump)/2 , rArray, nyArray))
        nyTestDummy.append(nyRespectR((riStart+rJump)/2 , rArray, nyArray))
    print("first part of nyTestDummy ", nyTestDummy[:20])
    print("last part of nyTestDummy ", nyTestDummy[-20:])
    # remove first and last element because they are known (sigma_r(r_center) and sigma_r(r_exterior))
    row.pop(0)
    col.pop(0)
    data.pop(0)
    row.pop()
    col.pop()
    data.pop()
    # add -sigma_r(r_exterior) to the constant vektor
    constant[-1] += -s_rOuter
    
    
    # second set of equations is assembled (one for each intervall of r) derived from eq. (13CA) with sigma_z=const.
    length = len(data)
    nyTestDummy2 = []
    for i, rJump in enumerate(rEnds):
        if i == 0:
            #print("windings, len(materialWidths), i, confinementDefined", windings, len(materialWidths), i, confinementDefined) ###########################################
            riStart = rCenter                                                                               
            constant.append(s_rCenter - (1+nyRespectR( (riStart+rJump)/2 , rArray, nyArray)) * calcIntegral(n=0, r=rJump, r_Start=riStart, r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)) 
        else:
            if confinementDefined == 1 and rJump == rEnds[-1]:
                jNew = jConfinement
                print("UND NUR EINMAL STEHT DAS HIER DA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            else:
                jNew = j
            riStart = rEnds[i - 1]
            constant.append(- (1+nyRespectR( (riStart+rJump)/2 , rArray, nyArray)) * calcIntegral(n=0, r=rJump, r_Start=riStart, r_Exterior=rExterior, r_Center=rCenter, j=jNew, b_za=b_za, b_zi=b_zi, b_0=b_0))
        # one extra row is added if the confinement is defined
        row.extend([windings*len(materialWidths) + i + confinementDefined,  windings*len(materialWidths) + i + confinementDefined,  windings*len(materialWidths) + i + confinementDefined,  windings*len(materialWidths) + i + confinementDefined])
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
    print(len(constant))
    # add -sigma_r(r_exterior) to the constant vektor
    constant[-1] += -s_rOuter
    #constant[-1] += 5000
    #print("Achtung, hier")
    #print(constant[-1])
    #print(constant[-2])
    #print(len(constant))
    
    # third set of equations is assembled (one for each material transition) derived from u(r_i-1End) = u_(r_iStart) with sigma_z=const.
    length = len(data)
    rEndsWithoutLastEnd = rEnds[:] # [:] needed to copy the list (otherwise a reference to the original list would be assigned)
    rEndsWithoutLastEnd.pop()
    print("len(rEnds) len(rEndsWithoutLastEnd)", len(rEnds), len(rEndsWithoutLastEnd))
    for i, rJump in enumerate(rEndsWithoutLastEnd):       # iterating over the material transitions ##########################
        if i == 0:
            riStart = rCenter
        else:
            riStart = rEnds[i - 1]
        # two extra rows were edded if the confinement is defined
        row.extend([2 * windings*len(materialWidths) + i + 2*confinementDefined, 2 * windings*len(materialWidths) + i + 2*confinementDefined, 2 * windings*len(materialWidths) + i + 2*confinementDefined])
        icol = 3*i + 1
        col.extend([icol, icol +1, icol +2])
        data.extend([1/nyRespectR( (rJump+riStart)/2, rArray, EArray),
                     -nyRespectR( (rJump+riStart)/2 , rArray, nyArray) / nyRespectR( (rJump+riStart)/2, rArray, EArray)
                    + nyRespectR( (rEnds[i+1]+rJump)/2 , rArray, nyArray) / nyRespectR( (rEnds[i+1]+rJump)/2, rArray, EArray), 
                    -1/nyRespectR( (rEnds[i+1]+rJump)/2, rArray, EArray)])###################################################
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
    #A = csr_matrix((data, (row, col)), shape=(3*windings*len(materialWidths) -1, 3*windings*len(materialWidths) -1), dtype=float)
    # one extra layer for the confinement
    A = csr_matrix((data, (row, col)), shape=(3*windings*len(materialWidths) -1 + confinementDefined*3, 3*windings*len(materialWidths) -1 + confinementDefined*3), dtype=float)
    #print(A.toarray())
    if windings < 20:
        print("row:", row)
        print("column", col)
        print("data", data)
    print("length of row, column and data", len(data), len(col), len(row))
    print("max value in row, and column", max(row), max(col))
    print("A.shape: ", A.shape)
    print("type(A): ", type(A))
    #print(constant)
    print("len(constant): ", len(constant))
    print("type(constant): ", type(constant))
    
    # calculate all Streses at the material transitions
    #print("A in sparse: ", A[windings*len(materialWidths) - 5:])
    #print("A in dense: ", A.todense())
    print("type(A.todense()): ", type(A.todense()))
    print("determinant: ", det(A.todense()), "equal to zero: ", det(A.todense()) == 0)
    
    
    #nonzero_indices = A.nonzero()
    nonzero_values = A.data
    #print(len(data), len(nonzero_values))
    
    #print(constant)
    #plt.scatter(col, row, c=nonzero_values, cmap="YlGnBu",)
    plt.scatter(col, row, c=nonzero_values, cmap="brg",)
    #plt.scatter(col, row, c=nonzero_values, cmap="turbo",)
    plt.grid(True)
    plt.colorbar()
    plt.ylabel("Zeile")
    plt.gca().invert_yaxis()
    plt.xlabel("Spalte")
    plt.show()
    
    
    Ainv = inv(A)
    #Ainv = spsolve(A, identity) # identity is from scipy.sparse
    result = Ainv.dot(constant)
    #print(result)
        
    
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
    
    
    ###Berechne u_r
    u_rArray =   ( 1/nyRespectR(rCenter, rArray, EArray) * ( - nyRespectR(rCenter, rArray, nyArray) * s_rArray + s_phiArray) ) * rDomains[0]
    e_phiArray = ( 1/nyRespectR(rCenter, rArray, EArray) * ( - nyRespectR(rCenter, rArray, nyArray) * s_rArray + s_phiArray) )
    e_rArray = ( 1/nyRespectR(rCenter, rArray, EArray) * ( - nyRespectR(rCenter, rArray, nyArray) * s_phiArray + s_rArray) )
    #e_rArray = ( 1/nyRespectR(rCenter, rArray, EArray) * ( - nyRespectR(rCenter, rArray, nyArray) * s_rArray + s_phiArray) )
    #e_phiArray = 1/rDomains[0] * s_phiArray
    e_zArray = (-nyRespectR(rCenter, rArray, nyArray) / nyRespectR(rCenter, rArray, EArray) ) * (s_rArray + s_phiArray)
    
    for layer in range(len(rDomains) - 1):
        #rNew = np.linspace()######################################################
        layer += 1
        
        if confinementDefined == 1 and layer == (len(rDomains) - 2):
            print("DAS SOLLTE NUR EINMAL DA STEHEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            jNew = jConfinement
        else:
            jNew = j
        # print(result[3*layer], rEnds[layer - 1], len(rDomains[layer]))
        # print(len(- ((1 + nyRespectR( (rEnds[layer]-rEnds[layer-1]/2), rArray, nyArray)) / 2) * calcIntegral(0, rDomains[layer], rDomains[layer][0], rExterior, rCenter, j, b_za, b_zi, b_0)))
        # print(len(-  (1/rDomains[layer]**2) * (1 - (1 + nyRespectR( (rEnds[layer]-rEnds[layer-1]/2), rArray, nyArray))/2 ) * calcIntegral(2, rDomains[layer], rDomains[layer][0], rExterior, rCenter, j, b_za, b_zi, b_0)))
        # print(len(+  (result[3*layer - 1]  *  (1/2)  *  (1 + (rEnds[layer - 1]**2/r**2)))))
        # print(+  (result[3*layer - 1]  *  (1/2)  ))
        # print((1 + (rEnds[layer - 1])))
        # print(len(r**2))
        
        s_rArrayNew = (  (1/2) * result[3*layer] * (1 - rEnds[layer - 1]**2/rDomains[layer]**2) 
                    - ((1 + nyRespectR( (rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray)) / 2) * calcIntegral(n=0, r=rDomains[layer], r_Start=rDomains[layer][0], r_Exterior=rExterior, r_Center=rCenter, j=jNew, b_za=b_za, b_zi=b_zi, b_0=b_0)
                    -  (1/rDomains[layer]**2) * (1 - (1 + nyRespectR( (rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray)) /2 ) * calcIntegral(n=2, r=rDomains[layer], r_Start=rDomains[layer][0], r_Exterior=rExterior, r_Center=rCenter, j=jNew, b_za=b_za, b_zi=b_zi, b_0=b_0)
                    +  (result[3*layer - 1]  *  (1/2)  *  (1 + (rEnds[layer - 1]**2/rDomains[layer]**2))))
        # result[3*layer -1] = s_r,layer,Start -> last Value for layer = layers*len(materials) - 1 because range is subtracting one und one was manually added
        
        
        s_phiArrayNew = ( (result[3*layer] - nyRespectR( (rEnds[layer]+rEnds[layer-1])/2 , rArray, nyArray) * s_zBegin)  
                        +  nyRespectR( (rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray) * s_z  
                        -  (1 + nyRespectR( (rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray)) * calcIntegral(n=0, r=rDomains[layer], r_Start=rDomains[layer][0], r_Exterior=rExterior, r_Center=rCenter, j=jNew, b_za=b_za, b_zi=b_zi, b_0=b_0)
                        + result[3*layer - 1]
                        - s_rArrayNew )
        
        u_rArrayNew = ( 1/nyRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, EArray) * ( - nyRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray) * s_rArrayNew + s_phiArrayNew) ) * rDomains[layer]
    
        # e_rArrayNew = ( 1/nyRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, EArray) * ( - nyRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray) * s_rArrayNew + s_phiArrayNew) )
        # #e_phiArrayNew = 1/rDomains[0] * s_phiArrayNew###############################
        # e_phiArrayNew = 1/rDomains[layer] * s_phiArrayNew
        e_rArrayNew = ( 1/nyRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, EArray) * ( - nyRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray) * s_phiArrayNew + s_rArrayNew) )
        #e_phiArrayNew = 1/rDomains[0] * s_phiArrayNew###############################
        e_phiArrayNew = ( 1/nyRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, EArray) * ( - nyRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray) * s_rArrayNew + s_phiArrayNew) )
        e_zArrayNew = (-nyRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray) / nyRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, EArray) ) * (s_rArrayNew + s_phiArrayNew)
    
        
        s_rArray = np.append(s_rArray, s_rArrayNew)
        s_phiArray = np.append(s_phiArray, s_phiArrayNew)
        u_rArray = np.append(u_rArray, u_rArrayNew)
        e_rArray = np.append(e_rArray, e_rArrayNew)
        e_phiArray = np.append(e_phiArray, e_phiArrayNew)
        e_zArray = np.append(e_zArray, e_zArrayNew)
        #print(layer)
    
    
    #print(nyRespectR(r_i, rArray, nyArray))
    #print(r_i)
    return s_rArray, s_phiArray, u_rArray, e_rArray, e_phiArray, e_zArray
    
    
# def calcAnaliticalDisplacement(rDomains, rCenter, s_r, s_phi, windings, EArray, nyArray):
#     rArray = rDomains.flatten()
#     rEnds = [item[0] for item in rDomains]
#     rEnds.pop(0)
#     rEnds.append(rDomains[-1][-1])
    
#     u_rArray = ( 1/nyRespectR(rCenter, rArray, EArray) * ( - nyRespectR(rCenter, rArray, nyArray) * s_r + s_phi) ) * rDomains[0]
    
#     for winding in range(len(rDomains) - 1):
#         winding += 1
    
#         u_rArrayNew = ( 1/nyRespectR((rEnds[winding]+rEnds[winding-1])/2, rArray, EArray) * ( - nyRespectR((rEnds[winding]+rEnds[winding-1])/2, rArray, EArray) * s_r + s_phi) ) * rDomains[winding]
#         u_rArray = np.append(u_rArray, u_rArrayNew)
    
#     return [u_rArray]
    
#### Letzter Wert in allen Domains, also Werte doppelt    
# def calcDomains(rCenter, materialWidths, numberOfWindings, numberOfValuesR):
#     '''Calculates an np.array containing arrays with discrete values coresponding to the material domains; only the last domain contains the end value'''
#     thisSpot = rCenter
#     results = []
#     for i in range(numberOfWindings):
#         for width in materialWidths:
#             results.append(np.linspace(thisSpot, thisSpot + width, int(numberOfValuesR/(numberOfWindings * len(materialWidths))) ))
#             thisSpot += width
#     return np.array(results)

def calcDomains(rCenter, materialWidths, numberOfWindings, numberOfValuesR):
    '''Calculates a list containing np.arrays with discrete values coresponding to the material domains; only the last domain contains the end value'''
    thisSpot = rCenter
    results = []
    for i in range(numberOfWindings):
        for testForLastValue, width in enumerate(materialWidths):
            domain = np.linspace(thisSpot, thisSpot + width, int(numberOfValuesR/(numberOfWindings * len(materialWidths))) )
            if not (i == (numberOfWindings - 1) and testForLastValue == (len(materialWidths) - 1)):
                domain = domain[:-1]
            results.append(np.array(domain))
            thisSpot += width
    return results

def calcMaterialValuesWithDomains(rDomains, materialValues):
    materialNys = materialValues # "Ny" is used as a dummy name, meaning any material parameter can be calculated
    "returns an array of Materialvalues corresponding to all r Values in rDomains"
    rDomainsNew = np.append(rDomains[:-1], [rDomains[-1][:-1]], axis=0)
    ny = np.ones_like(rDomainsNew)
    for i, nyValue in enumerate(materialNys):
        ny[i::len(materialNys)] = ny[i::len(materialNys)] * nyValue
    print(len(ny))
    print(type(ny))
    print("material.shape(): ", ny.shape, "\nrDomains without last value: rDomainsNew.shape(): ", rDomainsNew.shape, "\ntotal Number of Values: ", rDomainsNew.shape[0] * rDomainsNew.shape[1])
    return np.append(ny.flatten(), materialNys[-1])

def nyRespectR(r, rArray, nyArray):
    """calculates the materialparameter as an array with respect to the given r"""
    return np.interp(r, rArray, nyArray)

# used to calculate the displacement and strain according to caldwell
def calcDisplacement(r, s_r, s_phi, E, ny):
    #print(type(E), type(ny), type(s_r), type(s_phi), type(r))
    #print(np.size(E), np.size(ny), np.size(s_r), np.size(s_phi), np.size(r))
    u_r = ( 1/E * (-ny * s_r + s_phi) ) * r
    #print("ccccccccc")
    e_r = 1/E * (s_r - ny * s_phi)
    #e_r = r
    #print("test2")
    return([u_r, e_r])


######################################################################################
###### Main

#############################r = np.linspace(r_i, r_a, 500000)

rDom = calcDomains(r_i, materialWidths, totalWindings, numberOfValues)

print(np.shape(rDom[0]), np.shape(rDom[1]), np.shape(rDom[2]), "last length:", np.shape(rDom[-1]))
print(rDom[0][-1],rDom[1][0],rDom[1][-1],rDom[2][0])
#r = np.linspace(r_i, rDom[-1][-1], len())
Ny = calcMaterialValuesWithDomains(rDom, materialNys)
E = calcMaterialValuesWithDomains(rDom, materialEs)


# adding the confinement
if t_confinement != -1:
    rLastValue = rDom[-1][-1]
    rDom = rDom[:-1] + [rDom[-1][:-1]] # in rDom ist nur im letzten Eintrag der letzte r-Wert. Da eine Domain hinzugefügt wird, muss der letzte Eintrag entfernt werden
    r_confinement = np.linspace(rLastValue, rLastValue + t_confinement , int(numberOfValues/(totalWindings * len(materialWidths))) )
    print(type(rDom))
    rDom.append(r_confinement)
    if totalWindings < 20: 
        print("rDom: ", rDom)
    Ny = np.append(Ny[:-1], np.ones_like(r_confinement) * ny_confinement, axis=0)
    E = np.append(E[:-1], np.ones_like(r_confinement) * E_confinement, axis=0)


rDomFlattened = np.concatenate(rDom, axis=None).ravel()


print("size of rDomFlattened: ", rDomFlattened.size)
if totalWindings< 20:
    print("rDomFlattened: ", rDomFlattened)
    print("length of elements in rDom: ", [len(x) for x in rDom])
print("size of Ny: ", Ny.size)
print("Ny: ", Ny)
print("size of E: ", E.size)
print("E: ", E)

s_r, s_phi, u_rAnalytical, e_rAnalytical, e_phiAnalytical, e_zAnalytical = calcAnaliticalSolution(rDomains=rDom, rCenter=r_i, rExterior=r_a, s_rCenter=0, s_rOuter=0, s_zBegin=0, windings=totalWindings, nyArray=Ny, EArray=E, materialWidths=materialWidths, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)

### Berechnen der Verschiebungen
u_rCaldwell =  calcDisplacement(np.concatenate(rDomFlattened, axis=None).ravel(), calcStresses(r=rDomFlattened, r_a=r_a, r_i=r_i, s_z0=s_z0, s_ri=s_ri, s_ra=s_ra, nu=Ny[0], b_za=b_za, b_zi=b_zi, b_0=b_0, j=j)[0],
                                calcStresses(r=rDomFlattened, r_a=r_a, r_i=r_i, s_z0=s_z0, s_ri=s_ri, s_ra=s_ra, nu=Ny[0], b_za=b_za, b_zi=b_zi, b_0=b_0, j=j)[1],
                                materialEs[0], materialNys[0])[0]

#u_rAnalytical = calcAnaliticalDisplacement(rDomains=rDom, rCenter=r_i, s_r=s_r, s_phi=s_phi, windings=totalWindings, EArray=E, nyArray=Ny)[0]


# E(r)
plt.subplot(anzahlRowPlots, anzahlColumnPlots, 1)
plt.plot(rDomFlattened , E, label="vorgegebenes E(r)")
plt.xlabel(f"Radius in m E{int(np.log10(mSkalierung))}")
plt.ylabel(f"E-Modul in N/(m E{int(np.log10(mSkalierung))})^2")

# ny(r)
plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 1)
plt.plot(rDomFlattened , Ny, label="vorgegebenes ny(r)")
plt.xlabel(f"Radius in m E{int(np.log10(mSkalierung))}")
plt.ylabel("ny in 1")

plt.subplot(anzahlRowPlots, anzahlColumnPlots, 2)
plt.grid(True)
plt.plot(rDomFlattened, s_r)
plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 2)
plt.grid(True)
plt.plot(rDomFlattened, s_phi)


plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 3)
#u_rFourierFunction = calcDisplacement(r, s_rSolBvpFourier, s_phiSolBvpFourier, fourierFunctionE_f(r), fourierFunctionNy_f(r))[0]
#plt.plot(r, u_rFourierFunction, color="orange", label="u_r berechnet über BVP mit Fourierreihe",)
plt.plot(rDomFlattened, u_rAnalytical, label="u_r analytisch berechnet")
#plt.plot(r, u_rCaldwell, "b", label="u_r berechnet über Caldwell")
plt.plot(rDomFlattened, u_rCaldwell, "b", label="u_r berechnet über Caldwell")
plt.xlabel(f"Radius in m E{int(np.log10(mSkalierung))}")
plt.ylabel(f"u_r in m E{int(np.log10(mSkalierung))}")
plt.legend()

plt.subplot(anzahlRowPlots, anzahlColumnPlots, 4)
plt.plot(rDomFlattened, e_rAnalytical, label="e_r analytisch berechnet")
#plt.plot(r, e_rCaldwell * mSkalierung**(-1), "b", label="e_r berechnet über Caldwell")
#plt.plot(r, e_rDerivative, label=f"e_r berechnet über BVP mit analytischer Ableitung\ne_ri={round(e_rDerivative[0], 6)}; e_ra={round(e_rDerivative[-1], 6)}")
plt.xlabel(f"Radius in m E{int(np.log10(mSkalierung))}")
plt.ylabel(f"e_r in [-] E{int(np.log10(mSkalierung))}")
plt.legend()

plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 4)
#plt.plot(r, 1/r * mSkalierung**(-2) * calcStresses(r=r * mSkalierung**(-1), r_a=r_a * mSkalierung**(-1), r_i=r_i * mSkalierung**(-1), s_z0=s_z0, s_ri=s_ri, s_ra=s_ra, nu=ny[0], b_za=b_za, b_zi=b_zi, b_0=b_0, j=j * mSkalierung**2)[1], label="e_phi nach Caldwell")
plt.plot(rDomFlattened, e_phiAnalytical, label="e_phi analytisch berechnet")
#plt.plot(r, 1/r * s_phiSolBvpDerivativeTanh * mSkalierung**(-1), "--", label="e_phi berechnet als BVP mit analytischer Ableitung")
plt.ylabel("e_phi in [-]")
plt.legend()

plt.subplot(anzahlRowPlots, anzahlColumnPlots, 3)
plt.plot(rDomFlattened, e_zAnalytical, label="e_z analytisch berechnet")
#plt.plot(r, e_zCaldwell * mSkalierung**(-1), "b", label="e_z berechnet über Caldwell")
#plt.plot(r, e_zDerivative, label=f"e_z berechnet über BVP mit analytischer Ableitung\ne_zi={round(e_zDerivative[0], 6)}; e_za={round(e_zDerivative[-1], 6)}")
plt.xlabel("Radius in m")
plt.ylabel("e_z in [-]")
plt.legend()

# Save plot
fig = plt.gcf()
fig.set_size_inches(20, 12)
currentTime = datetime.now()
pictureName = f"Bilder\Graphen{currentTime.year}-{currentTime.month:02d}-{currentTime.day:02d}-{currentTime.hour:02d}-{currentTime.minute:02d}"
plt.savefig(pictureName, dpi=500)

#plt.tight_layout()
print("Fertig")
plt.show()