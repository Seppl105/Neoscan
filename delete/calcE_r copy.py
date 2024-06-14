import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def calcMaterialTanhCoefficients(r_innen, r_aussen, t, materialWidths, materialValues, steigung):
    numberOfMaterials = len(materialWidths)
    '''Berechnet diskret die differenzierbare Funktion E(r)'''

    class WindungsdickenError(Exception):
        '''Definierter Fehler mit Exeption'''
        pass

    def WindungsdickePrüfen(materialWidths, Windung):
        '''Prüft, ob die Summe der Dicken der drei einzelnen Layer gleich der Windungsdicke ist.'''
        totalWidth = 0
        for width in materialWidths:
            totalWidth += width
        if not abs(totalWidth - Windung) < 1e-9:
            raise WindungsdickenError(
                f"Die Summe der Dicken der drei Layer ist nicht gleich der Windungsdicke t: {Windung} != {totalWidth}.")
        else:
            print("Abmaße einer Windung und der einzelnen Layer sind verträglich.")

    # Sprungstellen nach jeder Schichtdicke
    def calcSprungstellen(r_innen, materialWidths, anzahlWindungen):
        aktuelleSprungstelle = r_innen
        results = []
        for i in range(anzahlWindungen):
            for width in materialWidths:
                aktuelleSprungstelle += width
                results.append(aktuelleSprungstelle)

        return np.array(results)

    def calcA(e):
        '''Berechnet A'''
        e = np.asarray(e)
        a = 1 / 2 * (e[1:] - e[:-1])
        return a

    def calcC(b, x):
        '''Berechnet C'''
        c = -1 * b[:] * x[:]
        return c

    def calcD(a):
        '''Berechnet D'''
        d = a
        return d

    # def rectangle(r, start, t_width):
    #     return np.where((r >= start - 4 * t_width) & (r <= start + 8 * t_width), 1.0, 0)
    #     *rectangle(x[:, None], abs(c / b), t_width=min(material))

    def sumJumps(x, a, b, c, d):
        '''Summiert über alle tanh-Funktionen'''
        print()
        return np.sum(a * np.tanh(np.outer(b * x) + c) + d, axis=1)

    try:
        WindungsdickePrüfen(materialWidths, t)
    except WindungsdickenError as e:
        print(e)

    # Windungszahl berechnen
    windungen = lambda x, y, z: (x - y) / z
    anzahlWindungen = int(windungen(r_aussen, r_innen, t))  # i.d.R sind es 600/q Windungen
    print("Anzahl der Windungen: ", anzahlWindungen)

    # f(x) = A * tanh(B*x + C) + D
    A = np.ones(numberOfMaterials * anzahlWindungen - 1)  # A[i] = DeltaE/2
    B = np.ones(numberOfMaterials * anzahlWindungen - 1) * steigung  # B[i] = >500 (fast) frei wählbar; muss groß genug sein
    C = np.ones(numberOfMaterials * anzahlWindungen - 1)  # C[i] = -B * r_Sprung
    D = np.ones(numberOfMaterials * anzahlWindungen - 1)  # D[i] = E_(i) + A

    # Drei E-Moduli pro Windung
    E = np.ones(numberOfMaterials * anzahlWindungen)
    for i in range(numberOfMaterials):
        E[i::numberOfMaterials] = materialValues[i]

    sprungstellen = calcSprungstellen(r_innen, materialWidths, anzahlWindungen)
    sprungstellen = sprungstellen[:-1]  # eine Sprungstelle zu viel

    A = np.array(calcA(E))
    C = calcC(B, sprungstellen)
    D = calcD(A)
    # print(A.size, B.size, C.size, D.size)
    print(A)
    # Diskretisierung des Radius' für alle Windungen
    r = np.linspace(r_innen, (r_innen + anzahlWindungen * t), 360 * anzahlWindungen)
    # Berechnung des E-Moduls über den Radius diskretisiert
    E_o = materialValues[0] + sumJumps(r, A, B, C, D)
    E = E_o + materialValues[0] - E_o[0] # nochmal nachvollziehen !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return r, E_o, E, [A, B, C, D]



def calcMaterialValue(r, coefficients, E0):
    def sumJumps(x, a, b, c, d, axis):
        '''Summiert über alle tanh-Funktionen'''
        return np.sum(a * np.tanh(b * x + c) + d, axis=axis)

    # print(np.ma.isarray(r), r)
    if not isinstance(r, np.ndarray):
        E = E0 + sumJumps(r, coefficients[0], coefficients[1], coefficients[2], coefficients[3], 0)
    else:
        E = []
        for value in r:
            # print("e")
            E = np.append(E,E0 + sumJumps(value, coefficients[0], coefficients[1], coefficients[2], coefficients[3], 0))
    return np.array(E)


def calcE_r(r, coeff, materialValues, E):
    E_r = np.zeros_like(r)
    for i in range(len(coeff[0])):
        E_r += coeff[0][i] * np.tanh(coeff[1][i] * r + coeff[2][i]) + coeff[3][i]
        # print(np.tanh(coeff[1][i] * r + coeff[2][i]))
        #print(coeff[1][i], coeff[2][i])
        #print(coeff[3][i])
    E_r += materialValues[0]
    # E_r = E_r - E[0] + materialValues[0]
    return E_r



# ### So wird die Funktion bisher aufgerufen

# r_i = 430 # [mm] innerer Radius
# r_a = 646 # [m] äußerer Radius

# t = 0.36 # [mm] Dicke einer Wicklung (Conductor, Cowinding, Insulation)
# mt = np.array([0.12,0.23,0.01])
# mE = np.array([500,450,400])
# s = 500

# r, E_o, E, [A, B, C, D] = calcMaterialTanhCoefficients(r_i,r_a,t,mt,mE,s)
# koeffizieten = [A, B, C, D]


# # Beispielhafte Berechnung an der Stelle r = []
# r_value = np.array([499.3, 500, 550, 603.4])
# E_r = calcE_r(r_value, koeffizieten, mE, E_o)
# # E_r = mt[0] + calcE_r(r=r_value, coeff=koeffizieten) - calcE_r(r=r_i, coeff=koeffizieten)
# print(r_value)
# print(E_r)




# # r = sp.symbols('r')
# # tanh_sum = sum(A[i] * sp.tanh(B[i] * r + C[i]) + D[i] for i in range(len(A)))
# # def evaluate_tanh_sum(r_value):
# #     return tanh_sum.subs(r, r_value)
# # r_value = 500
# # result = evaluate_tanh_sum(r_value)
# # print(result)


# # und geplottet
# # Plot E(r) für alle Windungen

# plt.figure(figsize=(8, 6))
# plt.plot(r, E_o, label='E(r) in N/mm^2', color='b')
# # Lable und Titel hinzufügen
# plt.xlabel('r')
# plt.ylabel('E(r)')
# plt.title('E(r) durch Tangens Hyperbolicus genähert')
# plt.grid(True)
# plt.legend()
# plt.show()
