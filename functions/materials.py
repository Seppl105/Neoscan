import numpy as np

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


    # def WindungsdickePrüfen(Conductor, Cowinding, Insulation, Windung):
    #     '''Prüft, ob die Summe der Dicken der drei einzelnen Layer gleich der Windungsdicke ist.'''
    #     if not abs((Conductor + Cowinding + Insulation) - Windung) < 1e-9:
    #         raise WindungsdickenError(
    #             f"Die Summe der Dicken der drei Layer ist nicht gleich der Windungsdicke t: {Windung} != {Conductor + Cowinding + Insulation}.")
    #     else:
    #         print("Abmaße einer Windung und der einzelnen Layer sind verträglich.")

    # Sprungstellen nach jeder Schichtdicke
    def calcSprungstellen(r_innen, materialWidths, anzahlWindungen):
        aktuelleSprungstelle = r_innen
        results = []
        for i in range(anzahlWindungen):
            for width in materialWidths:
                aktuelleSprungstelle += width
                results.append(aktuelleSprungstelle)
        
        return np.array(results)

    # def calcSprungstellen(r_i, con, cow, ins, n):

    #     results = []
    #     for i in range(anzahlWindungen):
    #         r_i += con
    #         results.append(r_i)
    #         r_i += cow
    #         results.append(r_i)
    #         r_i += ins
    #         results.append(r_i)
    #     return np.array(results)
    
    def calcA(e):
        '''Berechnet A'''
        e = np.asarray(e)
        a = 1 / 2 * (e[1:] - e[:-1])
        return a

    def calcC(b, x):
        '''Berechnet C'''
        c = -1 * b[:] * x[:]
        return c
    
    def calcD(e, a):
        '''Berechnet D'''
        d = e + a
        return d

    def sumJumps(x, a, b, c, d):
        '''Summiert über alle tanh-Funktionen'''
        return np.sum(a * np.tanh(b * x[:, None] + c) + d, axis=1)


    try:
        WindungsdickePrüfen(materialWidths, t)
    except WindungsdickenError as e:
        print(e)

    # Windungszahl berechnen
    windungen = lambda x, y, z: (x - y) / z
    anzahlWindungen = int(windungen(r_aussen, r_innen, t))  # i.d.R sind es 600/q Windungen
    print("Anzahl der Windungen: ",  anzahlWindungen )


    # f(x) = A * tanh(B*x + C) + D
    A = np.ones(numberOfMaterials * anzahlWindungen - 1)  # A[i] = DeltaE/2
    B = np.ones(numberOfMaterials * anzahlWindungen - 1) * steigung  # B[i] = >500 (fast) frei wählbar; muss groß genug sein
    C = np.ones(numberOfMaterials * anzahlWindungen - 1)  # C[i] = -B * r_Sprung
    D = np.ones(numberOfMaterials * anzahlWindungen - 1)  # D[i] = E_(i) + A/2
    # A = np.ones(3 * anzahlWindungen - 1)  # A[i] = DeltaE/2
    # B = np.ones(3 * anzahlWindungen - 1) * steigung  # B[i] = >500 (fast) frei wählbar; muss groß genug sein
    # C = np.ones(3 * anzahlWindungen - 1)  # C[i] = -B * r_Sprung
    # D = np.ones(3 * anzahlWindungen - 1)  # D[i] = E_(i) + A/2
    
    # Drei E-Moduli pro Windung
    E = np.ones(numberOfMaterials * anzahlWindungen)
    for i in range(numberOfMaterials):
        E[i::numberOfMaterials] = materialValues[i]
    # E = np.ones(3 * anzahlWindungen)
    # E[0::3] = E_con
    # E[1::3] = E_cow
    # E[2::3] = E_ins
    
    sprungstellen = calcSprungstellen(r_innen, materialWidths, anzahlWindungen)
    #sprungstellen = calcSprungstellen(r_innen, t_con, t_cow, t_ins, anzahlWindungen)
    sprungstellen = sprungstellen[:-1]  # eine Sprungstelle zu viel
    
    
    A = np.array(calcA(E))
    C = calcC(B, sprungstellen)
    D = calcD(E[:-1], abs(A))
    #print(A.size, B.size, C.size, D.size)
    # Diskretisierung des Radius' für alle Windungen
    r = np.linspace(r_innen, (r_innen + anzahlWindungen * t), 360 * anzahlWindungen)
    # Berechnung des E-Moduls über den Radius diskretisiert
    E = materialValues[0] + sumJumps(r, A, B, C, D)
    #E = E + materialValues[0] - E[0] # nochmal nachvollziehen !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return r, E, [A, B, C, D]

def calcMaterialValue(r, coefficients, E0):
    def sumJumps(x, a, b, c, d, axis):
        '''Summiert über alle tanh-Funktionen'''
        return np.sum(a * np.tanh(b * x + c) + d, axis=axis)
    #print(np.ma.isarray(r), r)
    if not isinstance(r, np.ndarray):
        E = E0 + sumJumps(r, coefficients[0], coefficients[1], coefficients[2], coefficients[3], 0)
    else:
        E = []
        for value in r:
            #print("e")
            E = np.append(E, E0 + sumJumps(value, coefficients[0], coefficients[1], coefficients[2], coefficients[3], 0))
    return np.array(E)


### So wird die Funktion bisher aufgerufen

# r,E = calcE_r(r_i,r_a,t,t_con,t_cow,t_ins,E_con,E_cow,E_ins,steigung=1000)
# und geplottet
# # Plot E(r) für alle Windungen
# plt.figure(figsize=(8, 6))
# plt.plot(r, E, label='E(r) in N/mm^2', color='b')
# # Lable und Titel hinzufügen
# plt.xlabel('r')
# plt.ylabel('E(r)')
# plt.title('E(r) durch Tangens Hyperbolicus genähert')
# plt.grid(True)
# plt.legend()
# plt.show()