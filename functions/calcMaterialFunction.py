import numpy as np
import matplotlib.pyplot as plt

def calcMaterialTanhCoefficients(r_innen, r_aussen, t, materialWidths, materialValues, slope=1000, scale=720):
    '''Calculates the coefficients of the smooth function E(r) that consists of hyperbolic tangents.'''
    '''Calculates at equally distributed spots in r the according YOUNG's Modulus.'''
    ### f(x) = A * tanh(B*x + C) + D
    
    t = t # t one winding in m
    materialWidths = np.array(materialWidths) # material width in m
    materialValues = np.array(materialValues) # material Values in SI-Basis Units
    # slope changes the factor infront of the variable inside of the tanh (B == slope)
    # scale changes the number of grid points
    
    class WindingmesureError(Exception):
        '''Definierter Fehler mit Exeption'''
        pass

    def checkWindingmeasure(materialWidths, Winding):
        '''Checks whether the sum of the thicknesses of the three individual layers is equal to the winding thickness.'''
        totalWidth = 0
        for width in materialWidths:
            totalWidth += width
        if not abs(totalWidth - Winding) < 1e-9:
            raise WindingmesureError(f"Die Summe der Dicken der drei Layer ist nicht gleich der Windungsdicke t: {Winding} != {totalWidth}.")
        else:
            print("Abmaße einer Windung und der einzelnen Layer sind verträglich.")

    def calcJumpspot(r_innen, materialWidths, numberOfWindings):
        '''Calculates all radii (results) at which the jumps occur.'''
        thisSpot = r_innen
        results = []
        for i in range(numberOfWindings):
            for width in materialWidths:
                thisSpot += width
                results.append(thisSpot)
        return np.array(results)

    def calcA(e):
        '''Calculates A'''
        # e = np.asarray(e)
        a = 1 / 2 * (e[1:] - e[:-1])
        return a

    def calcC(b, x):
        '''Calculates C'''
        c = -1 * b[:] * x[:]
        return c

    def calcD(a):
        '''Calculates D'''
        d = a
        return d

    def sumJumps(x, a, b, c, d):
        '''Sums up all tanh-functions'''
        return np.sum(a * np.tanh(np.outer(x, b) + c) + d, axis=1)

    try:
        checkWindingmeasure(materialWidths, t)
    except WindingmesureError as e:
        print(e)

    numberOfMaterials = len(materialWidths)

    # Calculate Number of Windings
    windings = lambda x, y, z: (x - y) / z
    numberOfWindings = int(windings(r_aussen, r_innen, t))  # should be 600/q Windings with given measurements
    print("Anzahl der Windungen: ", numberOfWindings)


    A = np.ones(numberOfMaterials * numberOfWindings - 1)  # A[i] = DeltaE/2
    B = np.ones(numberOfMaterials * numberOfWindings - 1) * slope  # B[i] = >500 (almost) free to choose; must be high enough
    C = np.ones(numberOfMaterials * numberOfWindings - 1)  # C[i] = -B * r_jump
    D = np.ones(numberOfMaterials * numberOfWindings - 1)  # D[i] = E_(i) + A

    # Three YOUNG's Modulus for each winding
    E = np.ones(numberOfMaterials * numberOfWindings)
    for i in range(numberOfMaterials):
        E[i::numberOfMaterials] = materialValues[i]

    jumpspot = calcJumpspot(r_innen, materialWidths, numberOfWindings)
    jumpspot = jumpspot[:-1]  # last spot not of interest

    A = calcA(E)
    C = calcC(B, jumpspot)
    D = calcD(A)

    # # Discretisation for all radii of every winding
    # r = np.linspace(r_innen, (r_innen + numberOfWindings * t), scale * numberOfWindings)
    # # Calculates YOUNG's Modulus discretised over the radius
    # E_o = materialValues[0] + sumJumps(r, A, B, C, D)

    return [A, B, C, D]


def calcTanhValue(r, coeff, materialValues):
    '''Calculates at any spots in r (r provided by slove.bvp) the according YOUNG's Modulus.'''
    E_r = np.zeros_like(r)
    E_r = np.sum(coeff[0] * np.tanh(np.outer(r, coeff[1]) + coeff[2]) + coeff[3], axis=1)
    # for i in range(len(coeff[0])):
    #     E_r += coeff[0][i] * np.tanh(coeff[1][i] * r + coeff[2][i]) + coeff[3][i]
    E_r += materialValues[0]
    if len(E_r) == 1:
        E_r = E_r[0]
    return E_r

def calcTanhDerivativeValue(r, coeff):
    '''Calculate the derivative for any r (skalar or np.array) according to d[ A*tanh(B*r+C)+D ]/dr = A * B * (1 - tanh^2(B*r+C)) '''
    dE_r = np.zeros_like(r)
    dE_r = np.sum(coeff[0] * coeff[1] *   ( 1- (np.tanh(np.outer(r, coeff[1]) + coeff[2]))**2 )   , axis=1)
    #dE_r = np.sum(coeff[0] * coeff[1] / np.cosh(np.outer(r, coeff[1]) + coeff[2]) + coeff[3])**2 )   , axis=1))
    # for i in range(len(coeff[0])):
    #     E_r += coeff[0][i] * np.tanh(coeff[1][i] * r + coeff[2][i]) + coeff[3][i]
    if len(dE_r) == 1:
        dE_r = dE_r[0]
    return dE_r



# ### So wird die Funktion bisher aufgerufen

# r_i = 430 # [mm] innerer Radius
# r_a = 646 # [m] äußerer Radius

# t = 0.36 # [mm] Dicke einer Wicklung (Conductor, Cowinding, Insulation)
# mt = np.array([0.12,0.23,0.01])
# mE = np.array([500,450,400])
# s = 1000

# coefficients = calcMaterialTanhCoefficients(r_i,r_a,t,mt,mE,s, scale=900)



# # Beispielhafte Berechnung an der Stelle r = []
# r_value = 499.3
# E_r = calcTanhValue(r_value, coefficients, mE)
# print(r_value)
# print(E_r)


# # Plot E(r) für alle Windungen

# plt.figure(figsize=(8, 6))
# plt.plot(np.linspace(r_i, r_a, 500000), calcTanhValue(np.linspace(r_i, r_a, 500000), coefficients, mE), label='E(r) in N/mm^2', color='b')
# # Lable und Titel hinzufügen
# plt.xlabel('r')
# plt.ylabel('E(r)')
# plt.title('E(r) durch Tangens Hyperbolicus genähert')
# plt.grid(True)
# plt.legend()
# plt.show()



# def calcMaterialValue(r, coefficients, E0):
#     def sumJumps(x, a, b, c, d, axis):
#         '''Summiert über alle tanh-Funktionen'''
#         return np.sum(a * np.tanh(b * x + c) + d, axis=axis)
#
#     # print(np.ma.isarray(r), r)
#     if not isinstance(r, np.ndarray):
#         E = E0 + sumJumps(r, coefficients[0], coefficients[1], coefficients[2], coefficients[3], 0)
#     else:
#         E = []
#         for value in r:
#             # print("e")
#             E = np.append(E,E0 + sumJumps(value, coefficients[0], coefficients[1], coefficients[2], coefficients[3], 0))
#     return np.array(E)


# r = sp.symbols('r')
# tanh_sum = sum(A[i] * sp.tanh(B[i] * r + C[i]) + D[i] for i in range(len(A)))
# def evaluate_tanh_sum(r_value):
#     return tanh_sum.subs(r, r_value)
# r_value = 500
# result = evaluate_tanh_sum(r_value)
# print(result)
