#functions for use in drag_calc

#libraries
import math
import numpy as np
import scipy.io
from scipy import special
from numpy.linalg import matrix_power

#getareanrm function for draggin
def getareanrm():
    #data file with macromodel plate areas (m^2) and plate normal vectors in spacecraft body-fixed coordinates data refers to             #Starlink v1.0 and body-fixed coordinate system as shown in Starlink_geometry_20200331.jpg reflectivities and materials are           #guesses for now, satellite bus is assumed completely flat... need to revisit this at some point
    Au = 196.96655
    #SiO2 = 60.0843
    SiO2 = 20.03
    Al = 26.981538

    #arrays for returns
    A = np.zeros(6)
    plt_M = np.zeros(6)
    bwspec = np.zeros(6)
    bwdiff = np.zeros(6)
    nrm_sbf = np.zeros((3,6))


    #****SOLAR PANELS: Front panel -Y (Shark Fin mode) or -Z (Open Book mode, sada_pos=pi/2; default setting used here, to be changed     #outside of this routine)
    nrm_sbf[:,0] = [0, 0, -1] 
    A[0] = 8.05*2.77
    bwspec[0] = 0.05
    bwdiff[0] = 0.3
    plt_M[0] = SiO2

    #****SOLAR PANELS: Back panel +Y (Shark Fin mode) or +Z (Open Book mode, sada_pos=0; default setting used here, to be changed         #outside of this routine)
    nrm_sbf[:,1] = [0, 0, 1] 
    A[1] = 8.05*2.77
    bwspec[1] = 0.2
    bwdiff[1] = 0.4
    plt_M[1] = Al
    
    #****BUS: Zenith-facing -Z
    nrm_sbf[:,2] = [0, 0, 1] 
    A[2] = 2.82*1.27
    bwspec[2] = 0.2
    bwdiff[2] = 0.4
    plt_M[2] = Al 

    #***BUS: Nadir-facing +Z
    nrm_sbf[:,3] = [0, 0, -1] 
    A[3] = 2.82*1.27
    bwspec[3] = 0.2
    bwdiff[3] = 0.4
    plt_M[3] = Al    

    #****VISORS: Satellite-facing +Y
    nrm_sbf[:,4] = [0, 1, 0] 
    A[4] = .681+.426
    bwspec[4] = 0.2
    bwdiff[4] = 0.4
    plt_M[4] = Al  

    #****VISORS: Outward-facing -Y
    nrm_sbf[:,5] = [0, -1, 0] 
    A[5] = .681+.426
    bwspec[5] = 0.2
    bwdiff[5] = 0.4
    plt_M[5] = Al     

    for i in range(len(nrm_sbf)):
        nrm_sbf[:,i] = nrm_sbf[:,i] / np.linalg.norm(nrm_sbf[:,i])


    return  [nrm_sbf,A,plt_M,bwspec,bwdiff]


#ecef_llh for draggin
def ecef_llh(r_ecef,converge):
    #lat, lon in degrees
    #height in km
    a = 6378.1370; # km
    b = 6356.75231425; # km
    e = math.sqrt((a**2 - b**2)/a**2); # km

    x = r_ecef[:,0] # km, first value
    y = r_ecef[:,1] # km, second value
    z = r_ecef[:,2] # km, third value
    
    lon = math.atan2(y,x)
    p = math.sqrt(x**2 + y**2)
    lat = math.atan(z/(p*(1-e**2)))

    error = converge + 1

    while error > converge:

        N       = a**2/math.sqrt((a**2)*(math.cos(lat))**2 + (b**2)*(math.sin(lat))**2)
        height  = p/math.cos(lat) - N
        latn_1  = lat
        lat     = math.atan(z/(p*(1-(e**2)*N/(N + height))))
        error   = (abs(lat - latn_1)) #max taken out, might be needed later

    lat = lat*180/math.pi
    lon = lon*180/math.pi
    height = height

    return[lat,lon,height]


#attitude_model_rotmat for draggin
def attitude_model_rotmat(recef,sunVector_ecef,solar_lon_sbf):
    # %   Calculates the ECEF=>SBF rotation matrix based on
    # % recef (position in ECEF coordinates), sunVector_ECEF 
    # % (vector pointing at the sun in ECEF coordinates), and
    # % solar_lon_sbf (the longitude of the sun wrt the SBF coordinates
    # % in the SBF-x and -y plane
    # %   The method transforms this information into a rotation matrix by
    # % constructing a set of coordinates 't' where t1 is zenith facing, 
    # % t2, is perpendicular to t1 and the sunVector, and t3 completes the 
    # % right-handed coordinate frame. The t=>SBF and t=>ECEF rotations are
    # % constructed and then combined to give ECEF=>SBF

    tsbf = np.zeros((3,3))
    # tecef = np.zeros((3,3))

    tsbf[:,0] = [0,0,-1] #% t1 = zenith
    # tsbf[:,1] = np.cross(tsbf[:,1],[cosd(solar_lon_sbf(1)),sind(solar_lon_sbf(1)),0]) #% t2 = perpendicular to t1 and sunVector
    # tsbf[:,1] = np.cross(tsbf[:,0],[math.cos(math.radians(solar_lon_sbf(1))),math.sin(math.radians(solar_lon_sbf(1))),0])
    tsbf[:,1] = np.cross(tsbf[:,0],[math.cos(math.radians(solar_lon_sbf)),math.sin(math.radians(solar_lon_sbf)),0])
    tsbf[:,1] = tsbf[:,1]/math.sqrt(sum(tsbf[:,1]**2,0)) # normalize
    tsbf[:,2] = np.cross(tsbf[:,0],tsbf[:,1]) #% t3 = completes the right-hand coordinate

    tecef = np.zeros((3,3))
     
    tecef[:,0] = recef[:]/np.sqrt(sum(recef[:]**2,0)) #% t1 = zenith
    tecef[:,1] = np.cross(tecef[:,0],sunVector_ecef[0,:]) #% t2 = perpendicular to t1 and sunVector
    tecef[:,1] = tecef[:,1]/math.sqrt(sum(tecef[:,1]**2,1)) # normalize
    tecef[:,2] = np.cross(tecef[:,0],tecef[:,1]) #% t3 = completes the right-hand coordinate

    ecef2sbf = tsbf@(np.conj(tecef))
    
    return [ecef2sbf,tsbf,tecef]


#SentmanCD_MassSpecies function for draggin
def SentmanCD_MassSpecies(V,Ai,Normi,Ni,Ta,alpha):

    # % Estimate coefficient of Drag by variation of accomodation coefficient (alpha)
    # % INPUTS:
    # %   V is the Satellite Velocity wrt a co-rotating atmosphere
    # %   Ai is a vector containing plate areas (m^2)
    # %   Normi is a matrix containing Satellite Body Fixed (SBF) unit normal vectors (3xN matrix)
    # %   Ni is a vector of the number density of He, O, N2, O2, Ar, H, and N (cm^-3)
    # %   Ta is the temperature of the atmosphere
    # %   alpha is the accomodation coefficient (guessed as 0.93 by Bowman, 2007 (AIAA))
    # % OUTPUTS:
    # %   CD is the coefficient of drag for the satellite
    # %   Aref is the reference area of the satellite

    # % Universal Gas Constant
    R = 8314.47215 #% J/(K*kmol)
    # % Atomic Mass of He, O, N2, O2, AR, H, N
    Mass = [4.002602,15.9994,28.0134,31.9988,39.948,1.0079,14.0067] #Make into array
    Mass = np.array(Mass)
    # Mass = map(np.float(Mass))

    # % Mean Atomic Mass of Incident Atmosphere
    Ma_ratio = Ni*Mass/sum(Ni*Mass)
    # %Ma = sum(Ni.*Mass)/sum(Ni); % g/mol
    # % Most Probable thermal velocity of Incident Atmospheric Molecules
    vp = np.sqrt( 2*R*Ta/Mass )

    Tw = 300 #% Kelvin (Guessing)

    Vi    = math.sqrt( sum( V**2 ) )
    # Ti    = Mass@(matrix_power(Vi,2)/(3@R))
    # Ti    = np.multiply(Mass,(matrix_power(Vi,2)/(3*R)))
    Ti    = np.multiply(Mass,((Vi**2)/(3*R)))
    S1    = Vi/vp #% speed ratio: satellite speed to most probable thermal speed of ambient molecules
    Q     = 1 + 1/(2*S1**2)
    Vr = np.sqrt(2/3)*Vi*np.sqrt(1 + alpha*(Tw/Ti - 1))

    # % Pre-allocate memory
    PSIi = 0
    Aref  = 0

    for i in range(len(Ai)): #needs translation work

        gamma = np.dot(V,Normi[:,i])/Vi #% It is assumed that |Normi| = 1

        if gamma > 0 :
            Aref = Aref + Ai[i]*gamma
    
        P = np.exp(-(gamma**2)*S1**2) / S1 
        Z = 1 + special.erf(gamma*S1) 

	    # % Sum the product of CD_i*Aref_i = PSIi
	    # %     the appropriate coefficient of drag for the entire satellite is:
	    # %     C_{D} = sum( C_{D,i}*A_{ref,i} ) / sum( A_{ref,i} )
        PSIi = PSIi + Ai[i]*sum( Ma_ratio*( P/math.sqrt(math.pi) + gamma*Q*Z + (gamma*Vr/(2*Vi))*(gamma*math.sqrt(math.pi)*Z + P) ) )

    C_d = PSIi / Aref; #not sure if done, matrix division is making me unsure...

    return [C_d, Aref]

#main function, big bad boy, draggin
def draggin(recef,vecef,Tmsis,Dmsis,sun_deg,lat,lon,height):
    t_gps = np.array([0])
    # recef = np.zeros((1,3))
    # recef[0,0] = 6778137
    
    # vecef = np.array([0,0,7300]) #(m/2)
    #MSIS placeholder
    Nmsis = np.array([.5,.95,0,0,0,0,0]) #Number Densities (or ratio to total number density) of He, O, N2, O2, AR, H, N
    # Tmsis = 700 #temp at altitude (K)
    # Dmsis = 1e-12 #nuetral mass denisty (kg/m^3)
    alpha = .93

    #Vector pointing from satellite to sun (this should be obtainable from Orekit, and in ECEF coords)
    #but for now,
    # approx. solar vector at 1 AU sun-earth distance and around Dec. solstice, at midnight UT

    # sun_deg = -23.45
    sdeg1 = math.radians(sun_deg)
    sdeg_array = np.zeros((1,3))
    sdeg_array[0,0] = -math.cos(sdeg1)
    sdeg_array[0,2] = math.sin(sdeg1)
    # sdeg_array = np.array([-math.cos(sdeg1),0,math.sin(sdeg1)])

    sunVector_ecef = 149597870700*sdeg_array - recef

    [nrm_sbf,area,plt_M,bwspec,bwdiff] = getareanrm()

    scmass = 260 #Starlink v1.0 mass (kg) this is from wikipedia(!)
    att_model = scipy.io.loadmat('attitude_model.mat')
    # att_model.keys()
    mean_slonsbf = att_model['mean_slonsbf']

    mean_ssada = att_model['mean_ssada']

    latbin=[]
    sltbin=[]


    slt = (t_gps/3600 + lon/15) % 24 #approx. local time
    latedges = np.arange(-54, 54, 1)
    sltedges = np.arange(0,24,.25)


    lat_s = np.digitize(lat,latedges)
    sltbin = np.digitize(slt,sltedges)
    # slt_s = np.digitize(slt,sltedges)

    latbin.append(lat_s)

    ecef2sbf = np.zeros((3,3,len(t_gps)))
    vrel2sbf = np.zeros((3,len(t_gps)))
    vrel_sbf = np.zeros((3,len(t_gps)))
    Cd = np.zeros((1,len(t_gps)))
    Aref = np.zeros((1,len(t_gps)))

    # Loop thru example values
    ecef2sbf = np.zeros((3,3,len(t_gps)))
    sada_pos_interp = np.zeros((1,len(t_gps)))

    for j in range(len(t_gps)):
	    # % Lookup attitude model to calculate ECEF->SBF rotation matrix when quaternion interpolation is unreliable
	    [ecef2sbf[:,:,j-1],tsbf,tecef] = attitude_model_rotmat(recef,sunVector_ecef,mean_slonsbf[sltbin[j-1]-1,latbin[j-1]-1])
	    #perfect, 100% accurate to matlab

	    sada_pos_interp[j] = mean_ssada[sltbin[j]-1,latbin[j]-1]

	    # % Calculate velocity of satellite relative to atmosphere (ignoring winds but including co-rotation) in SBF coords
	    # boof = ecef2sbf[:,:,j-1]*vecef[:,j-1]
	    aa = ecef2sbf[:,:,j-1]
	    cc = vecef
	    # vrel_sbf[:,j-1] = ecef2sbf[:,:,j-1]*vecef
	    vrel_sbf[:,j-1] = np.matmul(ecef2sbf[:,:,j-1],vecef) #correct
	
	    Normi = np.zeros((3,6))
	    Normi[:,0] = ([0,-math.cos(sada_pos_interp[j-1]),-math.sin(sada_pos_interp[j-1])])
	    Normi[:,1] = ([0,math.cos(sada_pos_interp[j-1]),math.sin(sada_pos_interp[j-1])])
	    Normi[:,[2,3,4,5]] = nrm_sbf[:,[2,3,4,5]]
	    # Normi is right 


        # % Orient solar panel based on sada_pod angle (in radians)
	    # % (this won't work with parfor loops; instead modify the inputs to calc_srp below)
	    # %satprop.nrm_sbf (:,1) = [ 0; -cos(sada_pos_interp(j)); -sin(sada_pos_interp(j))]; % solar panels
	    # %satprop.nrm_sbf (:,2) = [ 0;  cos(sada_pos_interp(j));  sin(sada_pos_interp(j))]; % back of solar panels


	    [Cd,Aref] = SentmanCD_MassSpecies(vrel_sbf[:,j-1],area,Normi,Nmsis,Tmsis,alpha)

    drag_ecef = -0.5*(Cd*Aref/scmass)*Dmsis*np.sqrt(sum(vecef**2,1))*vecef #% drag acceleration in ECEF coordinates (m/s)
    drag_sbf = -0.5*(Cd*Aref/scmass)*Dmsis*np.sqrt(sum(vrel_sbf**2,1))*vrel_sbf #% drag acceleration in SBF coordinates (m/s)

    return [drag_ecef,drag_sbf]

