import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import *
from scipy.sparse.linalg import inv
from scipy.linalg import det
from functions.Hoopstress_Calc02Copy import calcStresses
from datetime import datetime # Nur für die Bennenung der Grafiken

roundToDigit = 9 # if a value is rounded it will be to the "roundToDigit" digits

#####   Defintion der Konstanten Werte
windingDivisor = 1       # Es wird für 600/windingDivisor Windungen gerechnet
numberOfValues = int(3000000 / windingDivisor) # Anzahl der Werte des r Vektors
totalWindings = int(600/windingDivisor)

#   Subplots
anzahlRowPlots = 2
anzahlColumnPlots = 4

#   Spannungsrandbedingungen
# bei Skalierung zu beachten: N/m^2 E-3 = kg m/(s^2 m^2) E-3 = kg/(s^2 m) E-3 = kg/(s^2 mm)
s_ra = 0    # [Pa] Spannung in r-Richtung außen
s_ri = 0    # [Pa] Spannung in r-Richtung innen 
s_z0 = 0    # [Pa] Spannung in z-Richtung

#   Materialparameter des Leiterbands
t_con = 0.12 * 10**(-3) # [m] Dicke des Bandleiters
t_cow = 0.23 * 10**(-3) # [m] Dicke des Cowindings
t_ins = 0.01 * 10**(-3) # [m] Dicke der Isolation
t = t_con + t_cow + t_ins             # [m] Dicke einer Wicklung (Conductor, Cowinding, Insulation)
materialWidths = [t_con, t_cow, t_ins]

E_con = 500 * 10**9#280 * 10**9 # [Pa] E-Modul des Conductors
E_cow = 450 * 10**9#300 * 10**9 # [Pa] E-Modul des Cowindings
E_ins = 400 * 10**9#200 * 10**9 # [Pa] E-Modul der Insulation
materialEs = [E_con, E_cow, E_ins]

ny_con = 0.35 # Possion's Ratio Conductor
ny_cow = 0.3  # Possion's Ratio Cowinding
ny_ins = 0.4  # Possion's Ratio Insulation
materialNys = [ny_con, ny_cow, ny_ins]

# Caldwell rechnet mit aus dem Leiterband gemittelten Materialparametern
ECaldwell = 0
nyCaldwell = 0
for i, width in enumerate(materialWidths):
    ECaldwell += width * materialEs[i]
    nyCaldwell += width * materialNys[i]
ECaldwell = ECaldwell / sum(materialWidths)
nyCaldwell = nyCaldwell / sum(materialWidths)


#   Materialparameter des Confinements
t_confinement = -1#0.0005       # [m] Dicke des Confinements (für "-1" wird ohne Confinement gerechnet)
E_confinement = E_con        # [Pa] E-Modul des Confinements
ny_confinement = ny_con      # Possions's Ratio des Confinements

#   Konstanten des Versuchsaufbau
r_i = 0.430               # [m] innerer Radius
r_a = 0.646               # [m] äußerer Radius ohne Confinement
if r_a != (r_i + t * ( ((r_a - r_i) / t) / windingDivisor)):
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n r_a und t Werte mit der Anzahl der Windungen ergeben nicht dev vorgegegebenen äußeren Radius\n vorgegebener Wert für r_a: ", r_a, "\n mit Dicken berechneter Wert ((r_i + t * ( ((r_a - r_i) / t) / windingDivisor)): ", (r_i + t * ( ((r_a - r_i) / t) / windingDivisor)), "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
r_a = (r_i + t * ( ((r_a - r_i) / t) / windingDivisor))


j = 114.2 * 10**6                        # [A/m^2] Stromdichte
b_za = -2.5                              # [T] magnetische Flussdichte außen
b_zi = 14                                # [T] magnetische Flussdichte innen
b_0 = b_za - (b_za-b_zi)/(r_a-r_i) * r_a # [T] absolutes Glied der Geradengleichung für B(r)




#####   Funktionen

def calcIntegral(n, r, r_Start, r_Exterior, r_Center, j, b_za, b_zi, b_0):
    # in Caldwell definiertes Integral
    return j * (  (b_za - b_zi)/(r_Exterior - r_Center) * (1 / (n+2) ) * (r**(n+2) - r_Start**(n+2)) 
                 +  b_0                     * (1/  (n+1) ) * (r**(n+1) - r_Start**(n+1)))
    

def calcAnaliticalSolution(rDomains, rCenter, rExterior, s_rCenter, s_rOuter, s_zBegin, windings, nyArray, EArray, materialWidths, j, b_za, b_zi, b_0):
    
    ### Berechne s_z
    s_z = s_zBegin  # from eq (4) ################################################################ ist nicht fertig eingebaut
    print(rDomains[-1][-1],)
    rArray = np.concatenate(rDomains, axis=None) # flatten the rArray (maybe a flat rArray would have been better in the first place)
    print(rArray[-1])
    print(rArray[-5:])
    print(type(rDomains))
    print("len(rArray) and len(rDomains[0]) and len(rDomains[-1])", len(rArray), len(rDomains[0]), len(rDomains[-1]))
    rEnds = [item[0] for item in rDomains]
    rEnds.pop(0)
    rEnds.append(rArray[-1]) #rEnds.append(rExterior)
    
    # if confinement is defined some parts of the calculation have to be adadpted
    # wheater the calculation has to be adapted is controled with the variable "confinementDefined"
    jConfinement = 0
    testVerschiebungConfinementNachInnen = 0 # Möglichkeit das Confinement als Innere Layer einzubringen
    confinementDefined = 0
    if len(rDomains) != windings*len(materialWidths):
        confinementDefined = 1
        print("confinement ist definiert mit rArray[-1]: ", rArray[-1], "und rExterior: ", rExterior)
        print("rEnds[-5:] ", rEnds[-5:])
    #print(rEnds)
    
    # lists to assemble the sparse matrix
    row = []
    col = []  # column
    data = []
    constant = [] # list of constants of the system of equations
    nyTestDummy = [] # is used to check i ny is assembled the right way
    
    
    # first set of equations is assembled (one for each intervall of r) derived from eq. (14CA) with sigma_z=const.
    cB = lambda r, r1 : 1/2 * (1 - r1**2/r**2)
    cA = lambda r, r1 : 1/2 * (1 + r1**2/r**2)
    c = lambda r, r1, ny, jC: - 1/2 * (1+ny) * calcIntegral(n=0, r=r, r_Start=r1, r_Exterior=rExterior, r_Center=rCenter, j=jC, b_za=b_za, b_zi=b_zi, b_0=b_0) - 1/r**2 * (1 - (1+ny)/2 ) * calcIntegral(n=2, r=r, r_Start=r1, r_Exterior=rExterior, r_Center=rCenter, j=jC, b_za=b_za, b_zi=b_zi, b_0=b_0)

    for i, rJump in enumerate(rEnds):
        if i == 0:
            riStart = rCenter
            constant.append(cA(rJump, riStart) * s_rCenter + c(rJump, riStart, materialRespectR( (riStart+rJump)/2 , rArray, nyArray), j)) 
        else:
            riStart = rEnds[i - 1]
            if confinementDefined == 1 and rJump == rEnds[-1 - testVerschiebungConfinementNachInnen]:
                jNew = jConfinement
                print("Erstes Set an Gleichungen aufgestellt und Confinement beachtet. rJump:", rJump)
            else:
                jNew = j
            constant.append(c(rJump, riStart, materialRespectR( (riStart+rJump)/2 , rArray, nyArray), jNew))
        row.extend([i, i, i])
        icol = 3*i - 1
        col.extend([icol, icol +1, icol +3])
        data.extend([-cA(rJump, riStart), -cB(rJump, riStart), 1])
        #print(materialRespectR((riStart+rJump)/2 , rArray, nyArray))
        nyTestDummy.append(materialRespectR((riStart+rJump)/2 , rArray, nyArray))
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
    constant[-1] = constant[-1] - s_rOuter
    
    
    # second set of equations is assembled (one for each intervall of r) derived from eq. (13CA) with sigma_z=const.
    length = len(data)
    nyTestDummy2 = []
    for i, rJump in enumerate(rEnds):
        if i == 0:
            #print("windings, len(materialWidths), i, confinementDefined", windings, len(materialWidths), i, confinementDefined) ###########################################
            riStart = rCenter                                                                               
            constant.append(s_rCenter - (1+materialRespectR( (riStart+rJump)/2 , rArray, nyArray)) * calcIntegral(n=0, r=rJump, r_Start=riStart, r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)) 
        else:
            if confinementDefined == 1 and rJump == rEnds[-1- testVerschiebungConfinementNachInnen]:
                jNew = jConfinement
                print("Zweites Set an Gleichungen aufgestellt und confinement beachtet. rJump:", rJump)
            else:
                jNew = j
            riStart = rEnds[i - 1]
            constant.append(- (1+materialRespectR( (riStart+rJump)/2 , rArray, nyArray)) * calcIntegral(n=0, r=rJump, r_Start=riStart, r_Exterior=rExterior, r_Center=rCenter, j=jNew, b_za=b_za, b_zi=b_zi, b_0=b_0))
        # if the confinement is defined, one extra row is added. Therfore, the row indices have to be adjusted by adding "confinementDefined"
        row.extend([windings*len(materialWidths) + i + confinementDefined,  windings*len(materialWidths) + i + confinementDefined,  windings*len(materialWidths) + i + confinementDefined,  windings*len(materialWidths) + i + confinementDefined])
        icol = 3*i - 1
        col.extend([icol, icol +1, icol +2, icol +3])
        data.extend([-1, -1, 1, 1])
        nyTestDummy2.append(materialRespectR((riStart+rJump)/2 , rArray, nyArray))
        
        if nyTestDummy[i] != materialRespectR((riStart+rJump)/2 , rArray, nyArray): # nur zum prüfen von ny
            print("Error !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("first part of nyTestDummy2", nyTestDummy2[:20])
    print("last part of nyTestDummy2", nyTestDummy2[-20:])
    # remove first and last element of the newly added values because they are know (sigma_r(r_center) and sigma_r(r_exterior))
    row.pop(length)    # not (length - 1) because the first element added after length was assigned shall be poped
    col.pop(length)
    data.pop(length)
    row.pop()
    col.pop()
    data.pop()
    # add -sigma_r(r_exterior) to the constant vektor
    constant[-1] = constant[-1] - s_rOuter
    
    
    # third set of equations is assembled (one for each material transition) derived from u(r_i-1End) = u_(r_iStart) with sigma_z=const.
    length = len(data)
    rEndsWithoutLastEnd = rEnds[:] # [:] needed to copy the list (otherwise a reference to the original list would be assigned)
    rEndsWithoutLastEnd.pop()
    for i, rJump in enumerate(rEndsWithoutLastEnd):       # iterating over the material transitions
        if i == 0:
            riStart = rCenter
        else:
            riStart = rEnds[i - 1]
        # If the confinement is defined, two extra rows were added. Therfore the row indices have to be adjusted by adding "2*confinementDefined"
        row.extend([2 * windings*len(materialWidths) + i + 2*confinementDefined, 2 * windings*len(materialWidths) + i + 2*confinementDefined, 2 * windings*len(materialWidths) + i + 2*confinementDefined])
        icol = 3*i + 1
        col.extend([icol, icol +1, icol +2])
        data.extend([1/materialRespectR( (rJump+riStart)/2, rArray, EArray),
                     -materialRespectR( (rJump+riStart)/2 , rArray, nyArray) / materialRespectR( (rJump+riStart)/2, rArray, EArray)
                    + materialRespectR( (rEnds[i+1]+rJump)/2 , rArray, nyArray) / materialRespectR( (rEnds[i+1]+rJump)/2, rArray, EArray), 
                    -1/materialRespectR( (rEnds[i+1]+rJump)/2, rArray, EArray)])###################################################

        if nyTestDummy[i] != materialRespectR((riStart+rJump)/2 , rArray, nyArray): #nur zum prüfen von ny
            print("Error !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    
    constant.extend([0] * (len(rEnds) - 1))
    
    #### calculate the sparse matrix
    # one extra layer and therfore 3 equations (/rows) if the confinement is defined
    A = csr_matrix((data, (row, col)), shape=(3*windings*len(materialWidths) -1 + confinementDefined*3, 3*windings*len(materialWidths) -1 + confinementDefined*3), dtype=float)
    #print(A.toarray())
    if windings < 20:
        print("constant", constant)
        print("row:", row)
        print("column", col)
        print("data", data)
    print("length of row, column and data", len(data), len(col), len(row))
    print("max value in row, and column", max(row), max(col))
    print("A.shape: ", A.shape)
    print("type(A): ", type(A))
    print("len(constant): ", len(constant))
    print("type(constant): ", type(constant))
    print("constant[windings*len(materialWidths) + confinementDefined - 1]", constant[windings*len(materialWidths) + confinementDefined - 1])
    print("constant[windings*len(materialWidths)*2 + confinementDefined*2 - 2]", constant[windings*len(materialWidths)*2 + confinementDefined*2 - 1])
    print("type(A.todense()): ", type(A.todense()))
    print("determinant: ", det(A.todense()), "equal to zero (nicht aussagekräftig): ", det(A.todense()) == 0)
    
    # plot the sparse matrix
    nonzero_values = A.data
    #print(len(data), len(nonzero_values))
    #plt.scatter(col, row, c=nonzero_values, cmap="YlGnBu",)
    plt.scatter(col, row, c=nonzero_values, cmap="brg",)
    #plt.scatter(col, row, c=nonzero_values, cmap="turbo",)
    plt.grid(True)
    plt.colorbar()
    plt.ylabel("Zeile")
    plt.gca().invert_yaxis()
    plt.xlabel("Spalte")
    plt.show()
   
    # plot the inverse matrix
    # Ainv = inv(A)
    # #plt.spy(Ainv.todense())
    # #plt.show()
    # print(type(Ainv))
    # nonzero_valuesAinv = Ainv.data
    # rowAinv, colAinv = Ainv.nonzero()
    # plt.scatter(colAinv, rowAinv, c=nonzero_valuesAinv, cmap="brg",)
    # plt.grid(True)
    # plt.colorbar()
    # plt.ylabel("Zeile")
    # plt.gca().invert_yaxis()
    # plt.xlabel("Spalte")
    # plt.show()
    
    
    result = np.linalg.solve(A.todense(), constant) # np.linalg.solve is faster and more accurate than inv(.) and .dot(Ainv)
    
    # print("np.linalg.cond(A)", np.linalg.cond(A.todense())) # condition number
    # Ainv = inv(A)
    # #Ainv = spsolve(A, identity) # identity is from scipy.sparse
    # result = Ainv.dot(constant)
    

    #toleranceRelativeAllclose = 10**(-8)
    #print("A und inv(inv(A)) sind mit relativer tolerance= ", toleranceTest, "  fast gleich:", np.allclose(Adense, Ainvinv, rtol=toleranceRelativeAllclose))
    

    ### Berechnen der Spannungen im innersten Torus
    # Berechne s_r
    s_rArray = (  (1/2) * result[0] * (1 - rCenter**2/rDomains[0]**2)                               
                - ((1 + materialRespectR(rCenter, rArray, nyArray)) / 2) * calcIntegral(n=0, r=rDomains[0], r_Start=rDomains[0][0], r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                -  (1/rDomains[0]**2) * (1 - (1 + materialRespectR(rCenter, rArray, nyArray)) /2 ) * calcIntegral(n=2, r=rDomains[0], r_Start=rDomains[0][0], r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                +  (s_rCenter  *  (1/2)  *  (1 + (rCenter**2/rDomains[0]**2))))
    
    
    # Berechne s_phi
    s_phiArray = ( (result[0] - materialRespectR(rCenter, rArray, nyArray) * s_zBegin)  
                 +  materialRespectR(rCenter, rArray, nyArray) * s_z  
                 -  (1 + materialRespectR(rCenter, rArray, nyArray)) * calcIntegral(n=0, r=rDomains[0], r_Start=rDomains[0][0], r_Exterior=rExterior, r_Center=rCenter, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                 + s_rCenter
                 - s_rArray )
    
    
    # Berechne u_r
    u_rArray =   ( 1/materialRespectR(rCenter, rArray, EArray) * ( - materialRespectR(rCenter, rArray, nyArray) * s_rArray + s_phiArray) ) * rDomains[0]
    e_phiArray = ( 1/materialRespectR(rCenter, rArray, EArray) * ( - materialRespectR(rCenter, rArray, nyArray) * s_rArray + s_phiArray) )
    e_rArray = ( 1/materialRespectR(rCenter, rArray, EArray) * ( - materialRespectR(rCenter, rArray, nyArray) * s_phiArray + s_rArray) )
    e_zArray = (-materialRespectR(rCenter, rArray, nyArray) / materialRespectR(rCenter, rArray, EArray) ) * (s_rArray + s_phiArray)
    

    resultSphiStart = result[0::3]
    resultSphiEnd = result[1::3]
    resultSrStart = result[2::3] # does not include the inner BC for s_rInnen
    print("resultSphiStart[-10:]: ", resultSphiStart[-10:])
    print("resultSphiEnd[-10:]: ", resultSphiEnd[-10:])
    print("resultSrStart[-10:]: ", resultSrStart[-10:])
    
    for layer in range(len(rDomains) - 1):
        layer += 1
        
        if confinementDefined == 1 and layer == (len(rDomains) - 1 - testVerschiebungConfinementNachInnen):
            print("Berechnen der Werte für rDomains[-1] (mit Confinement)")
            jNew = jConfinement
        else:
            jNew = j
        
        # s_rArrayNew = (  (1/2) * resultSphiStart[layer] * (1 - rEnds[layer - 1]**2/rDomains[layer]**2) 
        #             - ((1 + materialRespectR( (rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray)) / 2) * calcIntegral(n=0, r=rDomains[layer], r_Start=rDomains[layer][0], r_Exterior=rExterior, r_Center=rCenter, j=jNew, b_za=b_za, b_zi=b_zi, b_0=b_0)
        #             -  (1/rDomains[layer]**2) * (1 - (1 + materialRespectR( (rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray)) /2 ) * calcIntegral(n=2, r=rDomains[layer], r_Start=rDomains[layer][0], r_Exterior=rExterior, r_Center=rCenter, j=jNew, b_za=b_za, b_zi=b_zi, b_0=b_0)
        #             +  (resultSrStart[layer-1]  *  (1/2)  *  (1 + (rEnds[layer - 1]**2/rDomains[layer]**2))))
        s_rArrayNew = (  (1/2) * result[3*layer] * (1 - rEnds[layer - 1]**2/rDomains[layer]**2) 
                    - ((1 + materialRespectR( (rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray)) / 2) * calcIntegral(n=0, r=rDomains[layer], r_Start=rDomains[layer][0], r_Exterior=rExterior, r_Center=rCenter, j=jNew, b_za=b_za, b_zi=b_zi, b_0=b_0)
                    -  (1/rDomains[layer]**2) * (1 - (1 + materialRespectR( (rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray)) /2 ) * calcIntegral(n=2, r=rDomains[layer], r_Start=rDomains[layer][0], r_Exterior=rExterior, r_Center=rCenter, j=jNew, b_za=b_za, b_zi=b_zi, b_0=b_0)
                    +  (result[3*layer - 1]  *  (1/2)  *  (1 + (rEnds[layer - 1]**2/rDomains[layer]**2))))
        # result[3*layer -1] = s_r,layer,Start -> last Value for layer = layers*len(materials) - 1 because range is subtracting one und one was manually added
        
        s_phiArrayNew = ( (result[3*layer] - materialRespectR( (rEnds[layer]+rEnds[layer-1])/2 , rArray, nyArray) * s_zBegin)                          +  materialRespectR( (rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray) * s_z  
                        -  (1 + materialRespectR( (rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray)) * calcIntegral(n=0, r=rDomains[layer], r_Start=rDomains[layer][0], r_Exterior=rExterior, r_Center=rCenter, j=jNew, b_za=b_za, b_zi=b_zi, b_0=b_0)
                        + result[3*layer - 1]
                        - s_rArrayNew )

        
        u_rArrayNew = ( 1/materialRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, EArray) * ( - materialRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray) * s_rArrayNew + s_phiArrayNew) ) * rDomains[layer]
    

        e_rArrayNew = ( 1/materialRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, EArray) * ( - materialRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray) * s_phiArrayNew + s_rArrayNew) )
        #e_phiArrayNew = 1/rDomains[0] * s_phiArrayNew###############################
        e_phiArrayNew = ( 1/materialRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, EArray) * ( - materialRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray) * s_rArrayNew + s_phiArrayNew) )
        e_zArrayNew = (-materialRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, nyArray) / materialRespectR((rEnds[layer]+rEnds[layer-1])/2, rArray, EArray) ) * (s_rArrayNew + s_phiArrayNew)
    
        
        s_rArray = np.append(s_rArray, s_rArrayNew)
        s_phiArray = np.append(s_phiArray, s_phiArrayNew)
        u_rArray = np.append(u_rArray, u_rArrayNew)
        e_rArray = np.append(e_rArray, e_rArrayNew)
        e_phiArray = np.append(e_phiArray, e_phiArrayNew)
        e_zArray = np.append(e_zArray, e_zArrayNew)
    
    return s_rArray, s_phiArray, u_rArray, e_rArray, e_phiArray, e_zArray
    

def calcDomains(rCenter, materialWidths, numberOfWindings, numberOfValuesR):
    #Calculates a list containing np.arrays with discrete values coresponding to the material domains; only the last domain contains the end value
    # result = [ [r_1i] , [r_2i] , ...., [r_ni]] with i as a "free index"
    # to perform accurate additions the values have to be rounded (to avoid floition-point error)
    thisSpot = rCenter
    results = []
    for i in range(numberOfWindings):
        for testForLastValue, width in enumerate(materialWidths):
            domain = np.linspace(thisSpot, round(thisSpot + width, roundToDigit), int(numberOfValuesR/(numberOfWindings * len(materialWidths))) )
            if not (i == (numberOfWindings - 1) and testForLastValue == (len(materialWidths) - 1)):
                domain = domain[:-1]
            results.append(np.array(domain))
            thisSpot = round(thisSpot + width, roundToDigit)
    return results


def calcMaterialValuesWithDomains(rDomains, materialValues):
    # returns an np array of materialvalues witch the n-th element corresponding to the n-th r Values in rDomains
    materialNys = materialValues # "Ny" is used as a dummy name -> any material parameter can be calculated
    
    rDomainsNew = np.append(rDomains[:-1], [rDomains[-1][:-1]], axis=0)
    ny = np.ones_like(rDomainsNew)
    for i, nyValue in enumerate(materialNys):
        ny[i::len(materialNys)] = ny[i::len(materialNys)] * nyValue
    print(len(ny))
    print(type(ny))
    print("material.shape(): ", ny.shape, "\nrDomains without last value: rDomainsNew.shape(): ", rDomainsNew.shape, "\ntotal Number of Values: ", rDomainsNew.shape[0] * rDomainsNew.shape[1])
    return np.append(ny.flatten(), materialNys[-1])


def materialRespectR(r, rArray, materialArray):
    # interpolates the materialparameter from the given materialArray with respect to the given r
    return np.interp(r, rArray, materialArray)


    
def calcStressesCaldwell(r, r_a, r_i, s_z0, s_ra, s_ri, nu, b_za, b_zi, b_0, j):
    ### Berechne s_z
    s_z = s_z0  # from eq (4)
    
    ### Berechne s_phii
    s_phii = (  2 / (1 - (r_i**2/r_a**2))  ) * ( s_ra  +  ((1 + nu) / 2) * calcIntegral(n=0, r=r_a, r_Start=r_i, r_Exterior=r_a, r_Center=r_i, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                                                +  (1/r_a**2) * (1 - (1 + nu)/2 ) * calcIntegral(n=2, r=r_a, r_Start=r_i, r_Exterior=r_a, r_Center=r_i, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                                                -  (s_ri  *  (1/2)  *  (1 + (r_i**2/r_a**2))))
    #print("s_phii: ", s_phii)
    #print("Term 1: ", s_ra  +  ((1 + nu) / 2) * calcIntegral(0, r_a, r_a, r_i, j, b_za, b_zi, b_0), 
    #      "\nTerm 2: ", +  (1/r_a**2) * (1 - (1 + nu)/2 ) * calcIntegral(2,r_a, r_a, r_i, j, b_za, b_zi, b_0),
    #      "\nTerm 3: ",  -  (s_ri  *  (1/2)  *  (1 + (r_i**2/r_a**2))))
    
    ###Berechne s_r(r)
    s_rArray = (  (1/2) * s_phii * (1 - r_i**2/r**2) 
                - ((1 + nu) / 2) * calcIntegral(n=0, r=r, r_Start=r_i, r_Exterior=r_a, r_Center=r_i, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                -  (1/r**2) * (1 - (1 + nu)/2 ) * calcIntegral(n=2, r=r, r_Start=r_i, r_Exterior=r_a, r_Center=r_i, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                +  (s_ri  *  (1/2)  *  (1 + (r_i**2/r**2))))
    #print("s_rArray: ", s_rArray[0:5])
    
    ### Berechne s_phi
    s_phiArray = ( (s_phii - nu * s_z0)  
                 +  nu * s_z  
                 -  (nu + 1) * calcIntegral(n=0, r=r, r_Start=r_i, r_Exterior=r_a, r_Center=r_i, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)
                 + s_ri
                 - s_rArray )
    
    return[s_rArray, s_phiArray]


def calcDisplacementCaldwell(r, s_r, s_phi, E, ny):
    # used to calculate the displacement and strain according to caldwell
    u_r = ( 1/E * (-ny * s_r + s_phi) ) * r
    e_r = 1/E * (s_r - ny * s_phi)
    e_phi = 1/E * (-ny * s_r + s_phi)
    e_z = -ny/E * (s_r + s_phi)
    return([u_r, e_r, e_phi, e_z])


######################################################################################
##### Main

rDom = calcDomains(r_i, materialWidths, totalWindings, numberOfValues)

print(np.shape(rDom[0]), np.shape(rDom[1]), np.shape(rDom[2]), "last length:", np.shape(rDom[-1]))
print(rDom[0][-1],rDom[1][0],rDom[1][-1],rDom[2][0], rDom[-1][-1])

Ny = calcMaterialValuesWithDomains(rDom, materialNys)
E = calcMaterialValuesWithDomains(rDom, materialEs)

# adding the confinement
if t_confinement != -1:
    rLastValue = rDom[-1][-1]
    rDom = rDom[:-1] + [rDom[-1][:-1]] # in rDom ist nur im letzten Eintrag der letzte r-Wert. Da eine Domain hinzugefügt wird, muss der letzte Eintrag entfernt werden bevor die Domain fürs Confinement angehängt wird
    endValue = rLastValue + t_confinement
    r_confinement = np.linspace(rLastValue, endValue, int(numberOfValues/(totalWindings * len(materialWidths))) )
    print("!!!!!!!!!!!!!!",rLastValue, endValue, r_confinement[-1])
    print(type(rDom))
    rDom.append(r_confinement)
    if totalWindings < 20: 
        print("rDom: ", rDom)
    Ny = np.append(Ny[:-1], np.ones_like(r_confinement) * ny_confinement, axis=0)
    E = np.append(E[:-1], np.ones_like(r_confinement) * E_confinement, axis=0)


rDomFlattened = np.concatenate(rDom, axis=None)


print("size of rDomFlattened: ", rDomFlattened.size)
if totalWindings< 20:
    print("rDomFlattened: ", rDomFlattened)
    print("length of elements in rDom: ", [len(x) for x in rDom])
print("size of Ny: ", Ny.size)
print("Ny: ", Ny)
print("size of E: ", E.size)
print("E: ", E)

s_r, s_phi, u_rAnalytical, e_rAnalytical, e_phiAnalytical, e_zAnalytical = calcAnaliticalSolution(rDomains=rDom, rCenter=r_i, rExterior=r_a, s_rCenter=0, s_rOuter=0, s_zBegin=0, windings=totalWindings, nyArray=Ny, EArray=E, materialWidths=materialWidths, j=j, b_za=b_za, b_zi=b_zi, b_0=b_0)

### Berechnungen für Caldwell
# unterschiedlicher r-Vektor je nach dem ob das Confinement definiert wurde
if len(rDom) != totalWindings*len(materialWidths):
    # Confinement ist definiert, muss also aus dem r-Vektor zur Berechnung der Lösung von Caldwell entfernt werden
    rDomFlattenedCaldwell = np.concatenate(rDom[:-1], axis=None)
else:
    # Confinement ist nicht definiert, es kann also der selbe r-Vektor zur Berechnung der Lösung nach Caldwell genutzt werden
    rDomFlattenedCaldwell = rDomFlattened
    
s_zCaldwell = s_z0 # aus Impulsbilanz in z-Richtung
s_rCaldwell, s_phiCaldwell = calcStressesCaldwell(r=rDomFlattenedCaldwell, r_a=r_a, r_i=r_i, s_z0=s_z0, s_ri=s_ri, s_ra=s_ra, nu=nyCaldwell, b_za=b_za, b_zi=b_zi, b_0=b_0, j=j)
u_rCaldwell, e_rCaldwell, e_phiCaldwell, e_zCaldwell =  calcDisplacementCaldwell(r=rDomFlattenedCaldwell, s_r=s_rCaldwell, s_phi=s_phiCaldwell, E=ECaldwell, ny=nyCaldwell)

# E(r)
plt.subplot(anzahlRowPlots, anzahlColumnPlots, 1)
plt.plot(rDomFlattened , E, label="vorgegebenes E(r)")
plt.plot(rDomFlattenedCaldwell, [ECaldwell] * len(rDomFlattenedCaldwell), label="gemitteltes E für Caldwell")
plt.xlabel(f"Radius in m")
plt.ylabel(f"E-Modul in N/(m)^2")
plt.legend()

# ny(r)
plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 1)
plt.plot(rDomFlattened , Ny, label="vorgegebenes ny(r)")
plt.plot(rDomFlattenedCaldwell, [nyCaldwell] * len(rDomFlattenedCaldwell), label="gemitteltes ny für Caldwell")
plt.xlabel("Radius in m")
plt.ylabel("ny in 1")
plt.legend()

plt.subplot(anzahlRowPlots, anzahlColumnPlots, 2)
plt.grid(True)
plt.plot(rDomFlattened, s_r, label="s_r analytische Lösung")
plt.plot(rDomFlattenedCaldwell, s_rCaldwell, label="s_r nach Caldwell")
plt.xlabel("Radius in m")
plt.ylabel("s_r in Pa")
plt.legend()

plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 2)
plt.grid(True)
plt.plot(rDomFlattened, s_phi, label="s_phi analytische Lösung")
plt.plot(rDomFlattenedCaldwell, s_phiCaldwell, label="s_phi nach Caldwell")
plt.xlabel("Radius in m")
plt.ylabel("s_phi in Pa")
plt.legend()

plt.subplot(anzahlRowPlots, anzahlColumnPlots, 3)
plt.plot(rDomFlattened, e_rAnalytical, label="e_r analytische Lösung")
plt.plot(rDomFlattenedCaldwell, e_rCaldwell, label="e_r berechnet nach Caldwell")
#plt.plot(r, e_rDerivative, label=f"e_r berechnet über BVP mit analytischer Ableitung\ne_ri={round(e_rDerivative[0], 6)}; e_ra={round(e_rDerivative[-1], 6)}")
plt.xlabel(f"Radius in m")
plt.ylabel(f"e_r in [-]")
plt.legend()

plt.subplot(anzahlRowPlots, anzahlColumnPlots, 4)
#plt.plot(r, 1/r * mSkalierung**(-2) * calcStressesCaldwell(r=r * mSkalierung**(-1), r_a=r_a * mSkalierung**(-1), r_i=r_i * mSkalierung**(-1), s_z0=s_z0, s_ri=s_ri, s_ra=s_ra, nu=nyCaldwell, b_za=b_za, b_zi=b_zi, b_0=b_0, j=j * mSkalierung**2)[1], label="e_phi nach Caldwell")
plt.plot(rDomFlattened, e_phiAnalytical, label="e_phi analytische Lösung")
plt.plot(rDomFlattenedCaldwell, e_phiCaldwell, label="e_phi nach Caldwell")
#plt.plot(r, 1/r * s_phiSolBvpDerivativeTanh * mSkalierung**(-1), "--", label="e_phi berechnet als BVP mit analytischer Ableitung")
plt.ylabel("e_phi in [-]")
plt.legend()

plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 3)
#u_rFourierFunction = calcDisplacement(r, s_rSolBvpFourier, s_phiSolBvpFourier, fourierFunctionE_f(r), fourierFunctionNy_f(r))[0]
#plt.plot(r, u_rFourierFunction, color="orange", label="u_r berechnet über BVP mit Fourierreihe",)
plt.plot(rDomFlattened, u_rAnalytical, label="u_r analytisch berechnet")
#plt.plot(r, u_rCaldwell, "b", label="u_r berechnet über Caldwell")
plt.plot(rDomFlattenedCaldwell, u_rCaldwell, "b", label="u_r berechnet über Caldwell")
plt.xlabel(f"Radius in m")
plt.ylabel(f"u_r in m")
plt.legend()

plt.subplot(anzahlRowPlots, anzahlColumnPlots, anzahlColumnPlots + 4)
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