import numpy as np

def calcE_r(r_innen, r_aussen, t, t_con, t_cow, t_ins, E_con, E_cow, E_ins, steigung):
    '''Berechnet diskret die differenzierbare Funktion E(r)'''

    class WindungsdickenError(Exception):
        '''Definierter Fehler mit Exeption'''
        pass

    def WindungsdickePrüfen(Conductor, Cowinding, Insulation, Windung):
        '''Prüft, ob die Summe der Dicken der drei einzelnen Layer gleich der Windungsdicke ist.'''
        if not abs((Conductor + Cowinding + Insulation) - Windung) < 1e-9:
            raise WindungsdickenError(
                f"Die Summe der Dicken der drei Layer ist nicht gleich der Windungsdicke t: {Windung} != {Conductor + Cowinding + Insulation}.")
        else:
            print("Abmaße einer Windung und der einzelnen Layer sind verträglich.")

    try:
        WindungsdickePrüfen(t_con, t_cow, t_ins, t)
    except WindungsdickenError as e:
        print(e)

    # Windungszahl berechnen
    windungen = lambda x, y, z: (x - y) / z
    anzahlWindungen = int(windungen(r_aussen, r_innen, t))  # i.d.R sind es 600/q Windungen


    # f(x) = A * tanh(B*x + C) + D
    A = np.ones(3 * anzahlWindungen - 1)  # A[i] = DeltaE/2
    B = np.ones(3 * anzahlWindungen - 1) * steigung  # B[i] = >500 (fast) frei wählbar; muss groß genug sein
    C = np.ones(3 * anzahlWindungen - 1)  # C[i] = -B * r_Sprung
    D = np.ones(3 * anzahlWindungen - 1)  # D[i] = E_(i) + A/2

    # Drei E-Moduli pro Windung
    E = np.ones(3 * anzahlWindungen)
    E[0::3] = E_con
    E[1::3] = E_cow
    E[2::3] = E_ins

    # Sprungstellen nach jeder Schichtdicke
    def calcSprungstellen(r_i, con, cow, ins, n):

        results = []
        for i in range(anzahlWindungen):
            r_i += con
            results.append(r_i)
            r_i += cow
            results.append(r_i)
            r_i += ins
            results.append(r_i)
        return np.array(results)

    sprungstellen = calcSprungstellen(r_innen, t_con, t_cow, t_ins, anzahlWindungen)
    sprungstellen = sprungstellen[:-1]  # eine Sprungstelle zu viel

    def calcA(e):
        '''Berechnet A'''
        e = np.asarray(e)
        a = 1 / 2 * np.sign(e[1:] - e[:-1]) * abs(e[1:] - e[:-1])
        return a
    A = np.array(calcA(E))

    def calcC(b, x):
        '''Berechnet B'''
        c = -1 * b[:] * x[:]
        return c

    C = calcC(B, sprungstellen)
    def calcD(e, a):
        '''Berechnet D'''
        d = e + a
        return d

    D = calcD(E[:-1], abs(A))

    def sumJumps(x, a, b, c, d):
        '''Summiert über alle tanh-Funktionen'''
        return np.sum(a * np.tanh(b * x[:, None] + c) + d, axis=1)

    # Diskretisierung des Radius' für alle Windungen
    r = np.linspace(r_innen, (r_innen + anzahlWindungen * t), 360 * anzahlWindungen)
    # Berechnung des E-Moduls über den Radius diskretisiert
    E = E_con + sumJumps(r, A, B, C, D)
    E = E + E_con - E[0] # nochmal nachvollziehen !!!

    return r, E


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