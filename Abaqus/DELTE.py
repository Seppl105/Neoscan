import numpy as np
import matplotlib.pyplot as plt
import sympy as smp

# Eigenschaften der Spule
r_i = 430 # [mm] innerer Radius
r_a = 646 # [mm] äußerer Radius
t = 0.36 # [mm] Dicke einer Wicklung (Conductor, Cowinding, Insulation)
t_con = 0.12 # [mm] Dicke des Bandleiters
t_cow = 0.23 # [mm] Dicke des Cowindings
t_ins = 0.01 # [mm] Dicke der Isolation
E_con = 100 #* 10**9 # E-Modul Conductor
E_cow = 90 #* 10**9 # E-Modul Cowinding
E_ins = 50 #* 10**9 # E-Modul Insulation

# Anteil der drei Materialien pro Wicklung bestimmen
anteilConductor = t_con / t
anteilCowinding = t_cow / t
anteilInsulation = t_ins / t

class WindungsdickenError(Exception):
    '''Definierter Fehler mit Exeption'''
    pass

def WindungsdickePrüfen(Conductor, Cowinding, Insulation, Windung):
    '''Prüft, ob die Summe der Dicken der drei einzelnen Layer gleich der Windungsdicke ist.'''
    if not abs((Conductor + Cowinding + Insulation) - Windung) < 1e-9:
        raise WindungsdickenError(f"Die Summe der Dicken der drei Layer ist nicht gleich der Windungsdicke t: {Windung} != {Conductor + Cowinding + Insulation}.")
    else:
        print("Abmaße einer Windung und der einzelnen Layer sind verträglich.")
try:
    WindungsdickePrüfen(t_con, t_cow, t_ins, t)
except WindungsdickenError as e:
    print(e)

# Windungsanzahl berechnen
windungen = lambda x, y, z: (x - y) / z

anzahlWindungen = int(windungen(r_a, r_i, t) / 300)

# f(x) = A * tanh(B*x + C) + D
A = np.ones(3 * anzahlWindungen -1) # A[i] = DeltaE/2
B = np.ones(3 * anzahlWindungen -1) * 100 # B[i] = 200 (fast) frei wählbar; muss groß genug sein
# B[2::3] = 100
C = np.ones(3 * anzahlWindungen -1) # C[i] = -B * r_Sprung
D = np.ones(3 * anzahlWindungen -1) # D[i] = E_(i) + A/2

E = np.ones(3 * anzahlWindungen)
E[0::3] = E_con
E[1::3] = E_cow
E[2::3] = E_ins


def calcSprungstellen(r_innen,con,cow,ins,n):

    results = []
    for i in range(anzahlWindungen):
        r_innen += con
        results.append(r_innen)
        r_innen += cow
        results.append(r_innen)
        r_innen += ins
        results.append(r_innen)
    return np.array(results)

sprungstellen = calcSprungstellen(r_i,t_con,t_cow,t_ins, anzahlWindungen)
sprungstellen = sprungstellen[:-1]

def calcA(e):
    e = np.asarray(e)
    a = 1/2 * np.sign(e[1:]-e[:-1]) * abs(e[1:] - e[:-1])
    return a

def calcC(b,x):
    c = -1 * b[:] * x[:]
    return c

def calcD(e,a):
    d = e[:] + a[:] / 2
    return d


A = np.array(calcA(E))

C = calcC(B,sprungstellen)

D = calcD(E[:-1], A)
print(len(A),len(C),len(D))
print(A)
def sumJumps(x,a,b,c,d):
    '''Berechnet einen Sprung an der Stelle r'''
    return np.sum(a * np.tanh(b*x[:, None] + c) + d, axis=1)

# Diskretisierung des Radius für alle Windungen
r = np.linspace(r_i, (r_i+10*t), 360*anzahlWindungen)

E = sumJumps(r,A,B,C,D)
# E = E_con + sumJumps(r,A,B,C,D)

# def jump(r,a,b,c,d):
#     '''Berechnet einen Sprung an der Stelle r'''
#     x = smp.Symbol("x",real=True)
#     result = 0
#
#     for element in a * smp.tanh(b*x + c) + d:
#         result += element
#     print(type(result))
#
#     result_f = smp.lambdify(x,result)
#     return result_f(r)

# for value in r:





# # Diskretisierung des Radius für die erste Windung
# r_w1 = np.linspace(430, 430.36, 360*anzahlWindungen)
# # Näherung der Sprungfunktion für die erste Windung
# E_w1 = ( 100 - 5 * np.tanh(200*r_w1 - (200*430.12)) - 20 * np.tanh(200*r_w1 - (200*430.35)) + 25 * np.tanh(200*r_w1 - (200*430.36)) ) * 10**3 # [N/mm^2] E-Modul



# jumps = np.ones(3 * anzahlWindungen - 1) * (A * np.tanh(B*r + C)) + D


# Näherung der Sprungfunktion für die erste Windung
# E = 100 + jumps


#Plot E(r) für die erste Windung
plt.figure(figsize=(8, 6))
plt.plot(r, E, label='E(r) in N/mm^2', color='b')

# Add labels and a title
plt.xlabel('r')
plt.ylabel('E(r)')
plt.title('E(r) durch Tangens Hyperbolicus genähert')
plt.grid(True)
plt.legend()
plt.show()