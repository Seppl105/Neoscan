import numpy as np
import matplotlib.pyplot as plt

def calcMaterialTanhCoefficients(r_innen, r_aussen, t, materialWidths, materialValues, slope):
    '''Calculates the coefficients of the smooth function E(r) that consists of hyperbolic tangents.'''
    '''Calculates at equally distributed spots in r the according YOUNG's Modulus.'''
    ### f(x) = A * tanh(B*x + C) + D

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

    # Sprungstellen nach jeder Schichtdicke
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
    numberOfWindings = int(windings(r_aussen, r_innen, t))  # should be 600 Windings with given measurements
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

    # Discretisation for all radii of every winding
    r = np.linspace(r_innen, (r_innen + numberOfWindings * t), 720 * numberOfWindings)
    # Calculates YOUNG's Modulus discretised over the radius
    E_o = materialValues[0] + sumJumps(r, A, B, C, D)

    return r, E_o, [A, B, C, D]


def calcE_r(r, coeff, materialValues, E):
    '''Calculates at any spots in r (r provided by slove.bvp) the according YOUNG's Modulus.'''
    E_r = np.zeros_like(r)
    for i in range(len(coeff[0])):
        E_r += coeff[0][i] * np.tanh(coeff[1][i] * r + coeff[2][i]) + coeff[3][i]
    E_r += materialValues[0]
    return E_r


### So wird die Funktion bisher aufgerufen

r_i = 430 # [mm] innerer Radius
r_a = 646 # [m] äußerer Radius

t = 0.36 # [mm] Dicke einer Wicklung (Conductor, Cowinding, Insulation)
mt = np.array([0.12,0.23,0.01])
mE = np.array([500,450,400])
s = 1000

r, E_o, [A, B, C, D] = calcMaterialTanhCoefficients(r_i,r_a,t,mt,mE,s)
coefficients = [A, B, C, D]


# Beispielhafte Berechnung an der Stelle r = []
r_value = np.array([499.3, 500, 550, 603.4])
E_r = calcE_r(r_value, coefficients, mE, E_o)
print(r_value)
print(E_r)


# und geplottet
# Plot E(r) für alle Windungen

plt.figure(figsize=(8, 6))
plt.plot(r, E_o, label='E(r) in N/mm^2', color='b')
# Lable und Titel hinzufügen
plt.xlabel('r')
plt.ylabel('E(r)')
plt.title('E(r) durch Tangens Hyperbolicus genähert')
plt.grid(True)
plt.legend()
plt.show()



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
