import numpy as np




r_i = 430 # [mm] innerer Radius
r_a = 646 # [m] äußerer Radius

t = 0.36 # [mm] Dicke einer Wicklung (Conductor, Cowinding, Insulation)
t_con = 0.12
t_cow = 0.23
t_ins = 0.01

anteilConductor = t_con / t
anteilCowinding = t_cow / t
anteilInsulation = t_ins / t

print(anteilConductor,anteilCowinding,anteilInsulation)

E_con = 100 #* 10**9 # [N/m^2] E-Modul Conductor
E_cow = 90 #* 10**9 # [N/m^2] E-Modul Cowinding
E_ins = 50 #* 10**9 # [N/m^2] E-Modul Insulation

windungen = lambda x, y, z: (x - y) / z
diskretisierung = lambda n, m: n / m

anzahlWindungen = int(windungen(r_a,r_i,t)) # [-] 600 Windungen
anzahlWerte = 10 * int(diskretisierung(t, min(t_con, t_cow, t_ins))) # [-] 360 Werte pro Windung
print(anzahlWerte)
print(anteilConductor * anzahlWerte,anteilCowinding*anzahlWerte,anteilInsulation*anzahlWerte)
print(anzahlWindungen, anzahlWerte)

r = np.linspace(r_i,(r_i + t), anzahlWerte) # array mit diskreten Radien
E = np.ones(anzahlWerte)
E[:int(anteilConductor * anzahlWerte)] = E_con
E[int(anteilConductor * anzahlWerte) : int((anteilConductor + anteilCowinding) * anzahlWerte)] = E_cow
E[int((anteilConductor + anteilCowinding) * anzahlWerte) : int((anteilConductor + anteilCowinding + anteilInsulation) * anzahlWerte)] = E_ins


print(E)

#print(r)











# abstand = np.linspace(1,36, 360)
# E = np.ones(360)
# print(abstand)
#
# E_con = 100 #* 10**9 # [N/m^2] E-Modul Conductor
# E_cow = 90 #* 10**9 # [N/m^2] E-Modul Cowinding
# E_ins = 50 #* 10**9 # [N/m^2] E-Modul Insulation
#
# for i in abstand:
#    if i < (t_con * 10000):
#        E[int(i)] = E[int(i)] * E_con
#        print(i)
#    elif i < ((t_con + t_cow) * 10000):
#        E[int(i)] = E[int(i)] * E_cow
#        print(i)
#    elif i < ((t_con + t_cow + t_ins) * 10000):
#        E[int(i)] = E[int(i)] * E_cow
#        print(i)
#
# print(E)






np.ones(anzahlWerte)