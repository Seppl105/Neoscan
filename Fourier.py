# Fuktionswert an der Stelle r über eine Fourierreihe aus den Fourierkoeffizienten bestimmen
def inverseFourier(r, coefficients, frequencies):
    n = len(coefficients)
    result = np.real(coefficients[0]) / n # konstanter Anteil
    
    for k in range(1, n//2):  # aus symmetriegründen nur die erste Hälfte der Freuqenzen zu betrachten
        amplitude = 2 * np.abs(coefficients[k]) / n  # Amplitude berechnen
        phase = np.angle(coefficients[k])  # Phase berechnen
        result += amplitude * smp.cos(2 * 3.141592653589793 * frequencies[k] * r + phase) # Nur Realteil interessant also bleibt nur der cos(.)
        
    return result


# Fuktionswert an der Stelle r über eine Fourierreihe aus diskreten Inputwerten bestimmen
def FourierSeries(r, inputFunction):
    fftInputFunction = np.fft.fft(inputFunction) # Koeffizienten der fft berechnet mit der Funktion fft aus dem fft package
    return inverseFourier(r, fftInputFunction, np.fft.fftfreq(len(inputFunction)))

def dFourierSeries(r, inputFunction):
    fftInputFunction = np.fft.fft(inputFunction)
    coefficients = fftInputFunction
    frequencies = np.fft.fftfreq(len(inputFunction))
    
    n = len(coefficients)
    result = np.real(coefficients[0]) / n # konstanter Anteil
    
    for k in range(1, n//2):  # aus symmetriegründen nur die erste Hälfte der Freuqenzen zu betrachten
        amplitude = 2 * np.abs(coefficients[k]) / n  # Amplitude berechnen
        phase = np.angle(coefficients[k])  # Phase berechnen
        result += - amplitude * 2 * 3.141592653589793 * coefficients[k] * smp.sin(2 * 3.141592653589793 * frequencies[k] * r + phase) # wegen i als Vorfaktor durch die Ableitung ist für den Realteil nur -sin(.) interessant
        
    return result


# input Function for solbvp with E(r) and ny(r) as Fourrier Series
def funcFourier(r, y):
    dSigmadr = 1/r * y[1]   - calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   - 1/r * y[0],
    dSigmaPfi = -1/r * y[1]   + calcBFeld(r, r_a, r_i, b_za, b_zi, b_0) * j   + 1/r * y[0]   + fourierFunctionNy_f(r) * ds_z   + dfourierFunctionNy_f(r) * (y[0] + s_z)   + dfourierFunctionE_f(r) * 1/fourierFunctionE_f(r) * (-fourierFunctionNy_f(r) * y[0]  + y[1]  + fourierFunctionNy_f(r) * s_z)   - (1 + fourierFunctionNy_f(r)) * radialForce_f(r)
    
    return np.vstack((dSigmadr,dSigmaPfi))


## Berechnen der differenzierbaren Funktionen E(r) und Ny(r) und ihrer ersten Ableitungen nach r
#print("Berechne Fourrierreihen für E, ny und R und die entsprechenden Ableitungen")

# E(r)
x = smp.Symbol("x", real=True)
# fourierFunctionE = FourierSeries(x, E) # A function dependant on the "Symbol" x
# #print(FourierFunctionE)
# fourierFunctionE_f = smp.lambdify(x, fourierFunctionE) # Convert FourierFunctionEr to a numerical function for plotting
# #FourierFunctionEr = FourierSeries(np.linspace(0, len(r), len(r)*5), Er)


# # ny(r)
# fourierFunctionNy = FourierSeries(x, ny) # A function dependant on the "Symbol" x
# fourierFunctionNy_f = smp.lambdify(x, fourierFunctionNy) # Convert FourierFunctionEr to a numerical function for plotting
# #FourierFunctionNy = FourierSeries(np.linspace(0, len(r), len(r)*5), ny)


# # dE(r)/dr
# dfourierFunctionE = smp.diff(fourierFunctionE, x)
# dfourierFunctionE_f = smp.lambdify(x, dfourierFunctionE)
# dFormulaFourierFunctionE = dFourierSeries(x, E)
# dFormulaFourierFunctionE_f = smp.lambdify(x, dFormulaFourierFunctionE)

# # dny(r)/dr
# dfourierFunctionNy = smp.diff(fourierFunctionNy, x)
# dfourierFunctionNy_f = smp.lambdify(x, dfourierFunctionNy)
# dFormulaFourierFunctionNy = dFourierSeries(x, ny)
# dFormulaFourierFunctionNy_f = smp.lambdify(x, dFormulaFourierFunctionNy)

###   solving bvp with solve_bvp 
print("mit solvebvp")
##   mit Fourrierreiehn
print("für Fourrierreihen")

# Sguess = np.zeros((len(r), len(r)), dtype=float)
# yInitialGuess = np.zeros((2, r.size))

# solBvpFourier = solve_bvp(funcFourier, bc, r, yInitialGuess)

# s_rSolBvpFourier = solBvpFourier.sol(r)[0]
# s_phiSolBvpFourier = solBvpFourier.sol(r)[1]
