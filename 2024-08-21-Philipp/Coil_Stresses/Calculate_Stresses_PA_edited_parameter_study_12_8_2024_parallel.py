import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix
#from scipy.sparse.linalg import inv
#from scipy.linalg import det
#from functions.Hoopstress_Calc02Copy import calcStresses
from datetime import datetime # Nur für die Bennenung der Grafiken
#from multiprocessing import Pool
import multiprocessing


#####   Functions

def calcIntegral(n, r, r_Start, r_Exterior, r_Center, j, b_za, b_zi, b_0):
    # The integral as defined in the Caldwell paper
    return j * (  (b_za - b_zi)/(r_Exterior - r_Center) * (1 / (n+2) ) * (r**(n+2) - r_Start**(n+2)) 
                 +  b_0                     * (1/  (n+1) ) * (r**(n+1) - r_Start**(n+1)))
    

def calcAnaliticalSolution(rDomains, rCenter, rExterior, s_rCenter, s_rOuter, s_zBegin, windings, nyArray, EArray, materialWidths, j, b_za, b_zi,r_a,r_i):
    #The main calcuation of the physical stresses/strains/displacements within the coil

    b_0 = b_za - (b_za-b_zi)/(r_a-r_i) * r_a # Absolute term of the linear equation for B(r) [T]

    
    ### calcuate s_z
    s_z = s_zBegin  # from eq (4) ################################################################ ist nicht fertig eingebaut
    # print(rDomains[-1][-1],)
    rArray = np.concatenate(rDomains, axis=None) # flatten the rArray (maybe a flat rArray would have been better in the first place)
    
    rEnds = [item[0] for item in rDomains]
    rEnds.pop(0)
    rEnds.append(rArray[-1]) #rEnds.append(rExterior)
    
    # if confinement is defined some parts of the calculation have to be adadpted
    # wheater the calculation has to be adapted is controled with the variable "confinementDefined"
    jConfinement = 0
    test_confinement_within_coil = 0 # Case of Confinement as inner layer
    confinementDefined = 0
    if len(rDomains) != windings*len(materialWidths):
        confinementDefined = 1
    
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
            if confinementDefined == 1 and rJump == rEnds[-1 - test_confinement_within_coil]:
                jNew = jConfinement
                #print("Erstes Set an Gleichungen aufgestellt und Confinement beachtet. rJump:", rJump)
            else:
                jNew = j
            constant.append(c(rJump, riStart, materialRespectR( (riStart+rJump)/2 , rArray, nyArray), jNew))
        row.extend([i, i, i])
        icol = 3*i - 1
        col.extend([icol, icol +1, icol +3])
        data.extend([-cA(rJump, riStart), -cB(rJump, riStart), 1])
        #print(materialRespectR((riStart+rJump)/2 , rArray, nyArray))
        nyTestDummy.append(materialRespectR((riStart+rJump)/2 , rArray, nyArray))

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
            if confinementDefined == 1 and rJump == rEnds[-1- test_confinement_within_coil]:
                jNew = jConfinement
                #print("Zweites Set an Gleichungen aufgestellt und confinement beachtet. rJump:", rJump)
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
    
    result = np.linalg.solve(A.todense(), constant) # np.linalg.solve is faster and more accurate than inv(.) and .dot(Ainv)

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
    
    
    for layer in range(len(rDomains) - 1):
        layer += 1
        
        if confinementDefined == 1 and layer == (len(rDomains) - 1 - test_confinement_within_coil):
            #print("Berechnen der Werte für rDomains[-1] (mit Confinement)")
            jNew = jConfinement
        else:
            jNew = j
        
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
    

def calcDomains(rCenter, materialWidths, numberOfWindings, numberOfValuesR,roundToDigit):
    #Calculates a list containing np.arrays with discrete values coresponding to the material domains; only the last domain contains the end value
    # result = [ [r_1i] , [r_2i] , ...., [r_ni]] with i as a "free index"
    # to perform accurate additions the values have to be rounded (to avoid floition-point error)
    thisSpot = rCenter
    dom_out = []
    winding_type = []
    for i in range(numberOfWindings):
        for width_index, width in enumerate(materialWidths):
            domain = np.linspace(thisSpot, round(thisSpot + width, roundToDigit), int(numberOfValuesR/(numberOfWindings * len(materialWidths))) )
            if not (i == (numberOfWindings - 1) and width_index == (len(materialWidths) - 1)):
                domain = domain[:-1]
            dom_out.append(np.array(domain))
            winding_type.extend([width_index]*len(domain))
            thisSpot = round(thisSpot + width, roundToDigit)
    return dom_out, np.array(winding_type)


def calcMaterialValuesWithDomains(rDomains, materialValues):
    # returns an np array of materialvalues witch the n-th element corresponding to the n-th r Values in rDomains
    materialNys = materialValues # "Ny" is used as a dummy name -> any material parameter can be calculated
    
    rDomainsNew = np.append(rDomains[:-1], [rDomains[-1][:-1]], axis=0)
    ny = np.ones_like(rDomainsNew)
    for i, nyValue in enumerate(materialNys):
        ny[i::len(materialNys)] = ny[i::len(materialNys)] * nyValue
    # print(len(ny))
    # print(type(ny))
    # print("material.shape(): ", ny.shape, "\nrDomains without last value: rDomainsNew.shape(): ", rDomainsNew.shape, "\ntotal Number of Values: ", rDomainsNew.shape[0] * rDomainsNew.shape[1])
    return np.append(ny.flatten(), materialNys[-1])


def materialRespectR(r, rArray, materialArray):
    # interpolates the materialparameter from the given materialArray with respect to the given r
    return np.interp(r, rArray, materialArray)




def add_confinement(t_confinement,rDom,totalWindings,materialWidths,Ny,E,num_r_samples,ny_confinement,E_confinement):

    # adding the confinement
    if t_confinement != -1:
        rLastValue = rDom[-1][-1]
        rDom = rDom[:-1] + [rDom[-1][:-1]] # Only the last entry in rDom contains the last r value. Since a domain is being added, the last entry needs to be removed before the domain for confinement is appended
        endValue = rLastValue + t_confinement
        r_confinement = np.linspace(rLastValue, endValue, int(num_r_samples/(totalWindings * len(materialWidths))) )
        rDom.append(r_confinement)
        if totalWindings < 20: 
            print("rDom: ", rDom)
        Ny = np.append(Ny[:-1], np.ones_like(r_confinement) * ny_confinement, axis=0)
        E = np.append(E[:-1], np.ones_like(r_confinement) * E_confinement, axis=0)

    rDomFlattened = np.concatenate(rDom, axis=None)

    return Ny,E,rDom,rDomFlattened


def calc_avg_material_parameters(materialWidths,materialEs,materialNys):
    # Caldwell caclucates with average values over the differetn materials over one turn
    ECaldwell = 0
    nyCaldwell = 0
    for i, width in enumerate(materialWidths):
        ECaldwell += width * materialEs[i]
        nyCaldwell += width * materialNys[i]
    ECaldwell = ECaldwell / sum(materialWidths)
    nyCaldwell = nyCaldwell / sum(materialWidths)

    return ECaldwell, nyCaldwell



def bin_outer_coil_radious(r_i,r_a,full_winding_thickness):
    # Bin the outer diameter according to the winding thicknesses
    if r_a != (r_i + full_winding_thickness * ( ((r_a - r_i) / full_winding_thickness))):
        print("Winding thickness and number of winding do not match the specified coil width (defined by outer-inner radius)")
    r_a = (r_i + full_winding_thickness * ( ((r_a - r_i) / full_winding_thickness)))
    return r_a



def combine_materials(t_con, t_cow, t_ins, E_con, E_cow, E_ins, ny_con, ny_cow, ny_ins):
    
    full_winding_thickness = t_con + t_cow + t_ins# [m] Full thickness of one turn (Conductor, Cowinding, Insulation)
    materialEs = [E_con, E_cow, E_ins]
    materialWidths = [t_con, t_cow, t_ins]
    materialNys = [ny_con, ny_cow, ny_ins]
    
    return full_winding_thickness, materialEs, materialWidths, materialNys



def run_calculation(parameter_tupel):
    
    (totalWindings, samples_per_winding_section, roundToDigit,
        t_con, t_cow, t_ins, t_confinement, 
        E_con, E_cow, E_ins, E_confinement, 
        ny_con, ny_cow, ny_ins, ny_confinement, 
        r_i, r_a, 
        j, 
        b_za, b_zi, s_ra, 
        s_ri, s_z0)                                                = set_model_parameters()
    
    
    num_r_samples = int(totalWindings*samples_per_winding_section*3) # Number of r values
    
    E_cow = parameter_tupel[1]
    t_cow = parameter_tupel[0]
    
    
    (full_winding_thickness, 
     materialEs, 
     materialWidths, 
     materialNys) = combine_materials(t_con, 
                                      t_cow, 
                                      t_ins, 
                                      E_con, 
                                      E_cow, 
                                      E_ins, 
                                      ny_con, 
                                      ny_cow, 
                                      ny_ins)
                                      
    ECaldwell, nyCaldwell = calc_avg_material_parameters(materialWidths,materialEs,materialNys)
    
    r_a =  bin_outer_coil_radious(r_i, 
                                  r_a, 
                                  full_winding_thickness)
    
    rDom, winding_type = calcDomains(r_i, 
                                     materialWidths, 
                                     totalWindings, 
                                     num_r_samples, 
                                     roundToDigit)
    
    
    
    Ny = calcMaterialValuesWithDomains(rDom, materialNys)
    
    E = calcMaterialValuesWithDomains(rDom, materialEs)
    
    Ny,E,rDom,rDomFlattened =  add_confinement(t_confinement,rDom,totalWindings,materialWidths,Ny,E,num_r_samples,ny_confinement,E_confinement)
    
    (s_r, 
     s_phi, 
     u_rAnalytical, 
     e_rAnalytical, 
     e_phiAnalytical, 
     e_zAnalytical) = calcAnaliticalSolution(rDomains=rDom, 
                                                rCenter=r_i, 
                                                rExterior=r_a, 
                                                s_rCenter=0, 
                                                s_rOuter=0, 
                                                s_zBegin=0, 
                                                windings=totalWindings, 
                                                nyArray=Ny, 
                                                EArray=E, 
                                                materialWidths=materialWidths, 
                                                j=j, 
                                                b_za=b_za, 
                                                b_zi=b_zi,
                                                r_a=r_a,
                                                r_i=r_i)
                                             
    # print(f"Eval cowinding e-modulusus: {E_cow/(1e9):.3f} GPa")
    # print(f"Eval cowinding thickness: {t_cow*(1e6):.3f} micron")
    
    return np.max(s_phi[winding_type==0])


def set_model_parameters():

    
    totalWindings = 500

    samples_per_winding_section = 5

    #num_r_samples = 3000000 # Anzahl der Werte of r Vektors

    roundToDigit = 9 # if a value is rounded it will be to the "roundToDigit" digits


    # Material parameters of the windings
    t_con = 120 * 10**(-6) # [m] Thickness of the conductor
    t_cow = 230 * 10**(-6) # [m] Thickness of the cowindings
    t_ins =  10 * 10**(-6) # [m] Thickness of the insulation
    # t_confinement = 0.0005 # [m] Thickness of the confinements (use "-1" to calculate without confinement)
    t_confinement = -1 # [m] Thickness of the confinements (use "-1" to calculate without confinement)
    
    E_con = 200 * 10**9#280 * 10**9 # [Pa] E-modulus of conductors
    E_cow = 450 * 10**9#300 * 10**9 # [Pa] E-modulus of cowindings
    E_ins = 3 * 10**9#200 * 10**9 # [Pa] E-modulus of the insulation
    E_confinement = E_con  # [Pa] E-modulus of confinements
    
    ny_con = 0.35 # Poisson's ratio of conductor
    ny_cow = 0.3  # Poisson's ratio of cowinding
    ny_ins = 0.4  # Poisson's ratio of insulation
    ny_confinement = ny_con  # Poisson's ratio of confinements
    
    # Constants of the coil setup
    # Geometrical properties
    r_i = 0.457 # [m] Inner radius of coil
    r_a = 0.665 # [m] Outer radius of coil without confinement
    # Electromagnetic properties:
    j = 370 * 10**6 # [A/m^2] Current density
    # j = 114.2 * 10**6 # [A/m^2] Current density
    b_za = -4.7 # [T] Magnetic flux density outside
    b_zi = 16.5  # [T] Magnetic flux density inside
    
    # Stress boundary conditions
    # Note for scaling: N/m^2 E-3 = kg m/(s^2 m^2) E-3 = kg/(s^2 m) E-3 = kg/(s^2 mm)
    s_ra = 0    # [Pa] Stress in r-direction outside
    s_ri = 0    # [Pa] Stress in r-direction inside
    s_z0 = 0    # [Pa] Stress in z-direction


    return totalWindings,samples_per_winding_section, roundToDigit, t_con, t_cow, t_ins, t_confinement, E_con, E_cow, E_ins, E_confinement, ny_con, ny_cow, ny_ins, ny_confinement, r_i, r_a, j, b_za, b_zi, s_ra, s_ri, s_z0



#### Main Execution: ####
if __name__ == "__main__":
        
     
    # Study the influence between cowinding emodulusus, cowinding thickness and maximal hoopstress as parallized grid search

    cowinding_e_moduluss = np.arange(100,800,step=100)*(1e9)
    cowinding_thicknesses = np.arange(50,600,step=100)*(1e-6)
    thickness_grid,e_mod_grid = np.meshgrid(cowinding_thicknesses,cowinding_e_moduluss)
    hoopstress_grid = np.zeros_like(thickness_grid)
    
    study_vals =  list(zip(thickness_grid.ravel(), e_mod_grid.ravel()))
    
    num_cores = multiprocessing.cpu_count()
    
    with multiprocessing.Pool(num_cores) as p:
        hoopstress_results = p.map(run_calculation, study_vals)

    
    hoopstress_results = [float(x) for x in hoopstress_results]

    hoopstress_grid = np.array(hoopstress_results).reshape(thickness_grid.shape)


    plt.close("all")

    plt.figure()
    plt.pcolormesh(thickness_grid * 1e6, e_mod_grid * 1e-9, hoopstress_grid * 1e-6, shading='auto', cmap='viridis')
    plt.colorbar(label='Value')  # Add a colorbar for reference
    contour_levels = np.arange(50,1000,step=50)
    contour = plt.contour(thickness_grid * 1e6, e_mod_grid * 1e-9, hoopstress_grid * 1e-6, colors='white', levels=contour_levels)
    plt.clabel(contour, inline=True, fontsize=8, colors='black')  # Label the contour lines
    plt.xlabel('Thickness(cowinding)[Micron]')
    plt.ylabel('E-mod(cowinding)[GPa]')
    plt.title('Hoopstress[MPa]')
    plt.grid(True)

    # Show the plot
    plt.show()

    
    
    #plot_results(2,4,rDomFlattened,rDomFlattenedCaldwell,s_rCaldwell,s_phiCaldwell,e_rCaldwell,e_phiCaldwell,u_rCaldwell)
    
    # # Get the current date
    # current_date = datetime.now().strftime("%Y-%m-%d")
    
    # # Create a filename with the date
    # filename = f'hoop_stress_{current_date}.txt'
    
    # # Save the array to a text file with the date in the filename
    # np.savetxt(filename, hoopstress_grid, delimiter=' ', fmt='%.2e', header='Hoop Stress Data')
