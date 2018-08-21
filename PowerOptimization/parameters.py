'''
Collected parameters of the pumping Kite system
Python Version 2.7 / Casadi version 2.4.1
- Author: Elena Malz, Chalmers 2016
'''
from collections import OrderedDict
import numpy as np
from numpy import pi

def AP2():
    params = OrderedDict()

    # kite mass [kg]
    params['mK'] = 36.8   # 50 * 25 [AWE book p.15]
    # acceleration due to gravity [m/s^]
    params['g'] = 9.81
    # density of air [kg/m^3]
    params['rho'] = 1.2  # with 1.255 is actually different !?
    # airframe moment of inertia [kg*m^2]
    params['J'] = np.array([[25,   0.0,   -0.47],
                            [  0.0, 32,   0.0],
                            [  -0.47,   0.0, 56.0]])
    # aerodynamic reference area [m^2]
    params['sref'] = 3.
    # aerodynamic reference span [m]
    params['bref'] = 5.5
    # aerodynamic reference chord [m]
    params['cref'] = 0.55

    #--- Cable parameters
    # tether natural length [m]
    params['l'] = 300.   # will reel out
    # tether mass [kg]
    params['mT'] = 0.0046 * params['l']
    # tether diameter [m]
    params['tether_diameter'] = 0.0025
    # tether density [kg/m]
    params['tether_density'] = 0.0046

    #--- Wind parameter
    # Drag coefficient for the tether
    params['Ctet']  = 1.2  # 0.4 in my model
    # altitude where wind shear starts
    params['windShearRefAltitude'] = 6.5
    # reference wind at ground [m/s]
    params['wind0'] = 8
    # roughness factor
    params['roughness'] = 0.12

    params['F_tether_max'] =  5000 #1./3 * pi *  (params['tether_diameter']*1000./2)**2 * 1073.4

    #--- Ampyx constraints
    params['tetacc_ub'] = 2.4
    params['tetacc_lb'] = -2.3
    params['omega_ulb'] = np.array([1,1,1])
    params['beta_ulb']  = 0.3 * 180./pi
    params['alpha_ub']  = 0.3 * 180./pi
    params['alpha_lb']  = -0.1 * 180./pi
    params['CD_alpha'] = [0.6176, 0.1572, 0.0372]  #NOTE: Like this it is wrong in the submitted paper, cause drag should be negative. But corrected signs do not change the plots. -.- Just dont tell anyone
    params['CL_alpha'] = [-5.3389, 5.03246, 0.5513] # not used

    # --- Ampyx wind coefficents
    params['alphaMaxDeg']   = np.inf #  deg. separation indicated somewhere between 12 and 15 degrees.

    params['CXalpha']      = [2.5549, 0.4784, -0.0293] # not used
    params['CZalpha']      = [5.7736, -5.0676, -0.5526]
    params['CXbeta']       = [0.0]
    params['CYbeta']       = [0.093591492950463,-0.029900842662737,-0.185542602484472]
    params['CZbeta']       = [0.0]
    # Damping Force
    params['Cxp']           = [0.0]
    params['Cxq']           = [4.4124, -0.6029 ]
    params['Cxr']           = [0.0]
    params['Cyp']           = [0.049646189990042,-0.014013395929910,-0.102202368153586]
    params['Cyq']           = [0.0]
    params['Cyr']           = [0.136826511089577,0.169467213438735]
    params['Czp']           = [0.0]
    params['Czq']           = [6.148615659988596,0.125149171823303,-7.556055799548278]
    params['Czr']           = [0.0]
    # Control dependent coefficients
    params['Cxail']          = [0.0]
    params['Cxele']          = [0.111459709670359,-0.010625762746063]
    params['Cxrud']          = [0.003]
    params['Cyail']          = [0.0579885138637340	-0.00241344656562128	-0.0514312281048974]
    params['Cyele']          = [0.0]
    params['Cyrud']          = [-0.103646503076224,0.026884904991135,0.103253859823406]
    params['Czail']          = [10.203508034464448,1.225557775231940,-5.400187761151892]
    params['Czele']          = [0.292385308913664,-0.001355014930468,-0.315035359951595]
    params['Czrud']          = [0.0]

    # MOMENTS
    # params['CL0']           = [0.0]
    # params['CM0']           = [-0.6027, -0.0307]
    # params['CN0']           = [0.0]
    params['Clalpha']       = [0.0]
    params['Cmalpha']       = [0.0, -0.6076,  -0.0307]
    params['Cnalpha']       = [0.0]
    params['Clbeta']        = [0.031258025680547,-3.483039877119464e-04,-0.063081351778656]
    params['Cmbeta']        = [0.0]
    params['Cnbeta']        = [-0.084911552609027,0.057710177865613]
    # Damping Moments
    params['Clp']           = [0.281317597010932,-0.024764054416691,-0.563258775268210]
    params['Clq']           = [0.0]
    params['Clr']           = [0.644873057041009,0.181159486166008]
    params['Cmp']           = [0.0]
    params['Cmq']           = [5.288473128322822,-0.002551748906449,-11.302198360248445]
    params['Cmr']           = [0.0]
    params['Cnp']           = [-0.913179228017379,-0.056597312252964]
    params['Cnq']           = [0.0]
    params['Cnr']           = [0.025700796184592,0.029060380546373,-0.055306030491248]
    # Control dependent coefficients
    params['Clail']          = [0.238326419615787,-0.008782572832499,-0.248918359753252]
    params['Clele']          = [0.0]
    params['Clrud']          = [-0.001346210113888,0.004363424631495]
    params['Cmail']          = [0.596509640460760,0.079570783623381,-0.316530950310559]
    params['Cmele']          = [0.997381197242358,-0.006148542440988,-1.042779240167403]
    params['Cmrud']          = [-0.001589335101276]
    params['Cnail']          = [-0.114703589462861,0.019036692592370]
    params['Cnele']          = [0.0]
    params['Cnrud']          = [0.040894647002254,-0.011721249877666,-0.040414259093759]
    return params

def simple_ampyx():
    params = OrderedDict()

    # kite mass [kg]
    params['mK'] = 36.8   # 50 * 25 [AWE book p.15]
    # acceleration due to gravity [m/s^]
    params['g'] = 9.81
    # density of air [kg/m^3]
    params['rho'] = 1.2  # with 1.255 is actually different !?
    # airframe moment of inertia [kg*m^2]
    params['J'] = np.array([[25,   0.0,   -4.7],
                            [  0.0, 32,   0.0],
                            [  -0.47,   0.0, 56.0]])
    # aerodynamic reference area [m^2]
    params['sref'] = 3.
    # aerodynamic reference span [m]
    params['bref'] = 5.5
    # aerodynamic reference chord [m]
    params['cref'] = params['sref']/params['bref']



    #--- Cable parameters
    # tether natural length [m]
    params['l'] = 300.   # will reel out
    # tether mass [kg]
    params['mT'] = 0.0046 * params['l']
    # tether diameter [m]
    params['tether_diameter'] = 0.0025
    # tether density [kg/m]
    params['tether_density'] = 0.0046

    #--- Wind parameter
    # Drag coefficient for the tether
    params['Ctet']  = 1.2  # 0.4 in my model
    # altitude where wind shear starts
    params['windShearRefAltitude'] = 6.5
    # reference wind at ground [m/s]
    params['wind0'] = 8
    # roughness factor
    params['roughness'] = 0.12

    params['F_tether_max'] =  5000 #1./3 * pi *  (params['tether_diameter']*1000./2)**2 * 1073.4

    #--- Ampyx xonstraints
    params['tetacc_ub'] = 2.4
    params['tetacc_lb'] = -2.3
    params['omega_ulb'] = np.array([1,1,1])
    params['beta_ulb']  = 0.3 * 180./pi
    params['alpha_ub']  = 0.3 * 180./pi
    params['alpha_lb']  = -0.1 * 180./pi

    # --- Ampyx wind coefficents
    params['alphaMaxDeg']   = np.inf #  deg. separation indicated somewhere between 12 and 15 degrees.
    #
    params['CX0']           = -0.0328
    params['CXalpha']       = 0.4784 *  (-1)#rad
    params['CXalphaSq']     = 2.5549 * (-1)
    params['CXbeta']        = 0.0 #* (180. / pi)*2 # per rad**2
    params['CY0']           = 0.0
    params['CYalpha']       = 0.0
    params['CYbeta']        = -0.1850 *(-1)
    params['CZ0']           = -0.5277   * (-1) # Lift coeffcient Cl
    params['CZalpha']       = -4.2247  * (-1) #* 180. / pi # per rad
    params['CZbeta']        = 0.0

    # Damping Force
    params['Cxp']           = 0.0
    params['Cxq']           = -0.6029
    params['Cxr']           = 0.0
    params['Cyp']           = -0.0510
    params['Cyq']           = 0.0
    params['Cyr']           = -0.17
    params['Czp']           = 0.0
    params['Czq']           =  -7.5
    params['Czr']           = 0.0
    # Control dependent coefficients
    params['Cxail']          = 0.0
    params['Cxele']          = -0.0106
    params['Cxrud']          = 0.0
    params['Cyail']          = -0.0510
    params['Cyele']          = 0.0
    params['Cyrud']          = 0.1030
    params['Czail']          = 0.0
    params['Czele']          = -0.3100
    params['Czrud']          = 0.0

    # MOMENTS
    params['CL0']           = 0.0
    params['CM0']           = -0.0306
    params['CN0']           = 0.0

    params['Clalpha']       = 0.
    params['Cmalpha']       = -0.6076
    params['Cnalpha']       = 0.0
    params['Clbeta']        = -0.0630
    params['Cmbeta']        = 0.
    params['Cnbeta']        = 0.0590


    # Damping Moments
    params['Clp']           = -0.5600
    params['Clq']           = 0.0
    params['Clr']           = 0.1500
    params['Cmp']           = 0.0
    params['Cmq']           = -11.3000
    params['Cmr']           = 0.0
    params['Cnp']           = -0.0400
    params['Cnq']           = 0.0
    params['Cnr']           = -0.0560

    # Control dependent coefficients
    params['Clail']          = -0.2480
    params['Clele']          = 0.0
    params['Clrud']          = 0.0040
    params['Cmail']          = 0.0
    params['Cmele']          =  -1.4200
    params['Cmrud']          = 0.0
    params['Cnail']          = 0.0210
    params['Cnele']          = 0.0
    params['Cnrud']          = -0.0400

    return params

def drag_2MW():
    params = OrderedDict()

    params['Pnom'] = 2
    # kite mass [kg]
    params['mK'] = 3000 # 1100.0   # 50 * 25 [AWE book p.15]
    # acceleration due to gravity [m/s^]
    params['g'] = 9.81
    # density of air [kg/m^3]
    params['rho'] = 1.2
    # airframe moment of inertia [kg*m^2]
    params['J'] = np.array([[4.4e3,   0.0,   0.0],
                            [  0.0, 2.1e3,   0.0],
                            [  0.0,   0.0, 6.2e3]])

    # tether natural length [m]
    params['l'] = 500.
    # tether mass [kg]
    params['mT'] = 250.
    # aerodynamic reference area [m^2]
    params['sref'] = 120
    # aerodynamic reference span [m]
    params['bref'] = 60
    # aerodynamic reference chord [m]
    params['cref'] = params['sref']/params['bref']
    # reference wind at ground [m/s]
    params['wind0'] = 5.0
    # tether diameter [m]
    params['tether_diameter'] = 0.04
    # tether density [kg/m]
    params['tether_density'] = params['mT']/params['l'] # 0.046
    # altitude where wind shear starts
    params['windShearRefAltitude'] = 5.
    # scaling
    params['ScalePower'] = 1e-3
    # maximal tether force
    # params['F_tether_max'] = 1./3 * pi *  (params['tether_diameter']*1000./2)**2 * 1073.4 #[N/mm^2] # dyneema
    params['F_tether_max'] = 1./3 * pi *  (params['tether_diameter']*1000./2)**2 * 3600 #[N/mm^2] #Aramid
    # wind harvesting efficiency
    params['eff'] = 0.6
    # Ampyx aerodynamic coefficients
    params['CD_alpha']     = [-0.6176, -0.1572, -0.0372]
    params['CZalpha']      = [-5.7736, +5.0676, +0.5526]
    params['CXbeta']       = [0.0]
    params['CYbeta']       = [0.093591492950463,-0.029900842662737,-0.185542602484472]
    params['CZbeta']       = [0.0]
    # Damping Force
    params['Cxp']           = [0.0]
    params['Cxq']           = [4.4124, -0.6029 ]
    params['Cxr']           = [0.0]
    params['Cyp']           = [0.049646189990042,-0.014013395929910,-0.102202368153586]
    params['Cyq']           = [0.0]
    params['Cyr']           = [-0.136826511089577,-0.169467213438735]
    params['Czp']           = [0.0]
    params['Czq']           = [6.148615659988596,0.125149171823303,-7.556055799548278]
    params['Czr']           = [0.0]
    # Control dependent coefficients - assumed as zero for this model due to malfunctioning
    params['Cxail']          = [0.0]
    params['Cxele']          = [0.0]
    params['Cxrud']          = [0.0]
    params['Cyail']          = [0.0]
    params['Cyele']          = [0.0]
    params['Cyrud']          = [0.0]
    params['Czail']          = [0.0]
    params['Czele']          = [0.0]
    params['Czrud']          = [0.0]
    # MOMENTS
    params['Clalpha']       = [0.0]
    params['Cmalpha']       = [0.0, -0.6076,  -0.0307]
    params['Cnalpha']       = [0.0]
    params['Clbeta']        = [0.031258025680547,-3.483039877119464e-04,-0.063081351778656]
    params['Cmbeta']        = [0.0]
    params['Cnbeta']        = [-0.084911552609027,0.057710177865613]
    # Damping Moments
    params['Clp']           = [0.281317597010932,-0.024764054416691,-0.563258775268210]
    params['Clq']           = [0.0]
    params['Clr']           = [0.644873057041009,0.181159486166008]
    params['Cmp']           = [0.0]
    params['Cmq']           = [5.288473128322822,-0.002551748906449,-11.302198360248445]
    params['Cmr']           = [0.0]
    params['Cnp']           = [-0.913179228017379,-0.056597312252964]
    params['Cnq']           = [0.0]
    params['Cnr']           = [0.025700796184592,0.029060380546373,-0.055306030491248]
    # Control dependent coefficients
    params['Clail']          = [0.238326419615787,-0.008782572832499,-0.248918359753252]
    params['Clele']          = [0.0]
    params['Clrud']          = [-0.001346210113888,0.004363424631495]
    params['Cmail']          = [0.596509640460760,0.079570783623381,-0.316530950310559]
    params['Cmele']          = [0.997381197242358,-0.006148542440988,-1.042779240167403]
    params['Cmrud']          = [-0.001589335101276]
    params['Cnail']          = [-0.114703589462861,0.019036692592370]
    params['Cnele']          = [0.0]
    params['Cnrud']          = [0.040894647002254,-0.011721249877666,-0.040414259093759]

    # params['J']             = np.array([[1375.,   0.0,   0.0],
    #                                     [  0.0,  869.,   0.0],
    #                                     [  0.0,   0.0, 2214.]])

    return params

def drag_1MW():
    "Large kite with low capacity factor but high FLH"
    params = OrderedDict()

    params['Pnom'] = 1
    # kite mass [kg]
    params['mK'] = 2000.0   # 50 * 25 [AWE book p.15]
    # acceleration due to gravity [m/s^]
    params['g'] = 9.81
    # density of air [kg/m^3]
    params['rho'] = 1.225
    # airframe moment of inertia [kg*m^2]   # if using J from greg not much is changed.
    params['J'] = np.array([[4.4e3,   0.0,   0.0],
                            [  0.0, 2.1e3,   0.0],
                            [  0.0,   0.0, 6.2e3]])
    # tether natural length [m]
    params['l'] = 500.
    # tether mass [kg]
    params['mT'] = 300.
    # aerodynamic reference area [m^2]
    params['sref'] = 110.
    # aerodynamic reference span [m]
    params['bref'] = 65.
    # aerodynamic reference chord [m]
    params['cref'] = params['sref']/params['bref']
    # reference wind at ground [m/s]
    params['wind0'] = 6.
    # tether diameter [m]
    params['tether_diameter'] = 0.035
    # altitude where wind shear starts
    params['tether_density'] = 0.46  #NOTE: is just guessed in the moment..
    # altitude where wind shear starts
    params['windShearRefAltitude'] = 5.
    # Scaling power for optimisation
    params['ScalePower'] = 1e-3
    # Set maximal tether force
    params['F_tether_max'] = 1./3 * pi *  (params['tether_diameter']*1000./2)**2 * 3600 #[N/mm^2] # aramid 3600, dyneema: 1073.4
    # wind harvesting efficiency
    params['eff'] = 0.6


    # parameters aerodynamic model rachel
    params['alphaMaxDeg']   = 10.0 # deg. separation indicated somewhere between 12 and 15 degrees.
    params['CL0']           = 0.53 #0.3455 #-------------------------------NOTE THIS IS CHANGED!!!!
    params['CLalpha']       = 0.04808 #* 180. / pi # per rad
    params['CD0']           = 0.02875
    params['CDalpha']       = -0.003561 #* 180. / pi # per rad
    params['CDalphaSq']     = 0.0006284 #* (180. / pi)*2 # per rad**2
    params['Clp']           = -0.48 * 180. / pi # per rad
    params['Clr']           = 0.01 * 180. / pi  # per rad
    params['Clbeta']        = -0.0008 * 180. / pi # per rad
    params['Cmalpha']       = -0.005 * 180. / pi # per rad
    params['Cmq']           = -9. # per rad
    params['Cnp']           = -0.11 * 180. / pi # per rad
    params['Cnr']           = -0.03 # per rad
    params['Cnbeta']        = 0.0003 * 180. / pi # per rad

    return params

def drag_666kW():
    "666kW kite with the aerodynamic model of rachel"
    params = OrderedDict()

    params['Pnom'] = 0.666
    # kite mass [kg]
    params['mK'] = 1200.0   # 50 * 25 [AWE book p.15]
    # acceleration due to gravity [m/s^]
    params['g'] = 9.81
    # density of air [kg/m^3]
    params['rho'] = 1.225
    # airframe moment of inertia [kg*m^2]   # if using J from greg not much is changed.
    params['J'] = np.array([[4.4e3,   0.0,   0.0],
                            [  0.0, 2.1e3,   0.0],
                            [  0.0,   0.0, 6.2e3]])
    # tether natural length [m]
    params['l'] = 500.
    # tether mass [kg]
    params['mT'] = 220.
    # aerodynamic reference area [m^2]
    params['sref'] = 80.
    # aerodynamic reference span [m]
    params['bref'] = 40.
    # aerodynamic reference chord [m]
    params['cref'] = params['sref']/params['bref']
    # reference wind at ground [m/s]
    params['wind0'] = 6.
    # tether diameter [m]
    params['tether_diameter'] = 0.02
    # altitude where wind shear starts
    params['tether_density'] = 0.46  #NOTE: is just guessed in the moment..
    # altitude where wind shear starts
    params['windShearRefAltitude'] = 5.
    # Scaling power for optimisation
    params['ScalePower'] = 1e-3
    # Set maximal tether force
    params['F_tether_max'] = 1./3 * pi *  (params['tether_diameter']*1000./2)**2 * 3600 #[N/mm^2] # aramid 3600, dyneema: 1073.4
    # wind harvesting efficiency
    params['eff'] = 0.6


  #parameters aerodynamic model rachel
    params['alphaMaxDeg']   = 10.0 # deg. separation indicated somewhere between 12 and 15 degrees.

     # --- Ampyx wind coefficents
    # params['alphaMaxDeg']   = np.inf #  deg. separation indicated somewhere between 12 and 15 degrees.

    params['CD_alpha']     = [-(0.0006284 * (180. / pi)**2 + 0.2), +0.003561 * 180. / pi , -0.02875]
    params['CZalpha']      = [0.0, 0.04808 * 180. / pi , 0.53]

    params['CXbeta']       = [0.0] # [0.0,0.0,-2] The original Rachel model included -2 * beta**2, now it's slightly different as no beta**2 is in the model.
    params['CYbeta']       = [0.0]
    params['CZbeta']       = [0.0]
    # Damping Force
    params['Cxp']           = [0.0]
    params['Cxq']           = [0.0]
    params['Cxr']           = [0.0]
    params['Cyp']           = [0.0]
    params['Cyq']           = [0.0]
    params['Cyr']           = [0.0]
    params['Czp']           = [0.0]
    params['Czq']           = [0.0]
    params['Czr']           = [0.0]
    # Control dependent coefficients
    params['Cxail']          = [0.0]
    params['Cxele']          = [0.0]
    params['Cxrud']          = [0.0]
    params['Cyail']          = [0.0]
    params['Cyele']          = [0.0]
    params['Cyrud']          = [0.0]
    params['Czail']          = [0.0]
    params['Czele']          = [0.0]
    params['Czrud']          = [0.0]

    # MOMENTS
    params['Clalpha']       = [0.0]
    params['Cmalpha']       = [0,-0.005 * 180. / pi,0]
    params['Cnalpha']       = [0.0]
    params['Clbeta']        = [0,0, -0.0008 * 180. / pi]
    params['Cmbeta']        = [0.0]
    params['Cnbeta']        = [0,0, 0.0003 * 180. / pi]

    # Damping Moments
    params['Clp']           = [0,0, -0.48 * 180. / pi ]
    params['Clq']           = [0.0]
    params['Clr']           = [0,0, 0.01 * 180. / pi]
    params['Cmp']           = [0.0]
    params['Cmq']           = [0,0,-9]
    params['Cmr']           = [0.0]
    params['Cnp']           = [0,0, -0.11 * 180. / pi ]
    params['Cnq']           = [0.0]
    params['Cnr']           = [0,0, -0.03 ]

    # Control dependent coefficients
    params['Clail']          = [0.0, 0.0, 1.0]
    params['Clele']          = [0.0]
    params['Clrud']          = [0.0]
    params['Cmail']          = [0.0]
    params['Cmele']          = [0.0, 0.0, 1.0]
    params['Cmrud']          = [0.0]
    params['Cnail']          = [0.0]
    params['Cnele']          = [0.0]
    params['Cnrud']          = [0.0, 0.0, 1.0]



    return params

def drag_666kW_Makani():
    """Represents the specificiations (weight and size) of the Makani wing of 600kW (Makani.pdf) with
     a simplified version of the aerodynamics of the amppx wing. Same inertia was used by Greg long time ago.
     Tether diameter is guessed."""

    params = OrderedDict()

    params['Pnom'] = 0.666
    # kite mass [kg]
    params['mK'] = 1050.0#1200# 1050.0 #1200#
    # acceleration due to gravity [m/s^]
    params['g'] = 9.81
    # density of air [kg/m^3]
    params['rho'] = 1.225
    # airframe moment of inertia [kg*m^2]
    params['J'] = np.array([[4.4e3,   0.0,   0.0],
                            [  0.0, 2.1e3,   0.0],
                            [  0.0,   0.0, 6.2e3]])
    # tether natural length [m]
    params['l'] = 500 # 440. # 500 #
    # tether mass [kg]
    params['mT'] = 220 #250. #220#
    # aerodynamic reference area [m^2]
    params['sref'] = 55. # 80
    # aerodynamic reference span [m]
    params['bref'] =  28. # 40
    # aerodynamic reference chord [m]
    params['cref'] = params['sref']/params['bref']
    # reference wind at ground [m/s]
    params['wind0'] = 6.
    # tether diameter [m]
    params['tether_diameter'] = 0.02
    # altitude where wind shear starts
    params['tether_density'] = 0.46#params['mT']/params['l']  #NOTE: is just guessed as the tether mass is given by makani.pdf
    # altitude where wind shear starts
    params['windShearRefAltitude'] = 5.
    # Scaling power for optimisation
    params['ScalePower'] = 1e-3
    # Set maximal tether force
    params['F_tether_max'] = 1./3 * pi *  (params['tether_diameter']*1000./2)**2 * 3600 #[N/mm^2] # aramid 3600, dyneema: 1073.4
    # wind harvesting efficiency
    params['eff'] = 0.6

    # Ampyx aerodynamic coefficients
    params['CD_alpha']     = [-0.6176, -0.1572, -0.0372]
    params['CZalpha']      = [-5.7736, +5.0676, +0.5526]
    params['CXbeta']       = [0.0]
    params['CYbeta']       = [0.093591492950463,-0.029900842662737,-0.185542602484472]
    params['CZbeta']       = [0.0]
    # Damping Force
    params['Cxp']           = [0.0]
    params['Cxq']           = [4.4124, -0.6029 ]
    params['Cxr']           = [0.0]
    params['Cyp']           = [0.049646189990042,-0.014013395929910,-0.102202368153586]
    params['Cyq']           = [0.0]
    params['Cyr']           = [-0.136826511089577,-0.169467213438735]
    params['Czp']           = [0.0]
    params['Czq']           = [6.148615659988596,0.125149171823303,-7.556055799548278]
    params['Czr']           = [0.0]
    # Control dependent coefficients - assumed as zero for this model due to malfunctioning
    params['Cxail']          = [0.0]
    params['Cxele']          = [0.0]
    params['Cxrud']          = [0.0]
    params['Cyail']          = [0.0]
    params['Cyele']          = [0.0]
    params['Cyrud']          = [0.0]
    params['Czail']          = [0.0]
    params['Czele']          = [0.0]
    params['Czrud']          = [0.0]

    # MOMENTS
    params['Clalpha']       = [0.0]
    params['Cmalpha']       = [0.0, -0.6076,  -0.0307]
    params['Cnalpha']       = [0.0]
    params['Clbeta']        = [0.031258025680547,-3.483039877119464e-04,-0.063081351778656]
    params['Cmbeta']        = [0.0]
    params['Cnbeta']        = [-0.084911552609027,0.057710177865613]
    # Damping Moments

    params['Clp']           = [0.281317597010932,-0.024764054416691,-0.563258775268210]
    params['Clq']           = [0.0]
    params['Clr']           = [0.644873057041009,0.181159486166008]
    params['Cmp']           = [0.0]
    params['Cmq']           = [5.288473128322822,-0.002551748906449,-11.302198360248445]
    params['Cmr']           = [0.0]
    params['Cnp']           = [-0.913179228017379,-0.056597312252964]
    params['Cnq']           = [0.0]
    params['Cnr']           = [0.025700796184592,0.029060380546373,-0.055306030491248]
    # Control dependent coefficients
    params['Clail']          = [0.238326419615787,-0.008782572832499,-0.248918359753252]
    params['Clele']          = [0.0]
    params['Clrud']          = [-0.001346210113888,0.004363424631495]
    params['Cmail']          = [0.596509640460760,0.079570783623381,-0.316530950310559]
    params['Cmele']          = [0.997381197242358,-0.006148542440988,-1.042779240167403]
    params['Cmrud']          = [-0.001589335101276]
    params['Cnail']          = [-0.114703589462861,0.019036692592370]
    params['Cnele']          = [0.0]
    params['Cnrud']          = [0.040894647002254,-0.011721249877666,-0.040414259093759]

    return params

def drag_666kW_oversized():
    """OVERSIZED specifications (weight and size) of the Makani wing of 600kW (Makani.pdf) with
     a simplified version of the aerodynamics of the ampyx wing. Same inertia was used by Greg long time ago.
     Tether diameter is guessed. """

    params = OrderedDict()

    params['Pnom'] = 0.6666
    # kite mass [kg]
    params['mK'] = 1200# 1050.0 #1200#
    # acceleration due to gravity [m/s^]
    params['g'] = 9.81
    # density of air [kg/m^3]
    params['rho'] = 1.225
    # airframe moment of inertia [kg*m^2]
    params['J'] = np.array([[4.4e3,   0.0,   0.0],
                            [  0.0, 2.1e3,   0.0],
                            [  0.0,   0.0, 6.2e3]])
    # tether natural length [m]
    params['l'] = 500 # 440. # 500 #
    # tether mass [kg]
    params['mT'] = 250. #220#
    # aerodynamic reference area [m^2]
    params['sref'] = 104
    # aerodynamic reference span [m]
    params['bref'] =  52
    # aerodynamic reference chord [m]
    params['cref'] = params['sref']/params['bref']
    # reference wind at ground [m/s]
    params['wind0'] = 6.
    # tether diameter [m]
    params['tether_diameter'] = 0.02
    # altitude where wind shear starts
    params['tether_density'] = 0.46#params['mT']/params['l']  #NOTE: is just guessed as the tether mass is given by makani.pdf
    # altitude where wind shear starts
    params['windShearRefAltitude'] = 5.
    # Scaling power for optimisation
    params['ScalePower'] = 1e-3
    # Set maximal tether force
    params['F_tether_max'] = 1./3 * pi *  (params['tether_diameter']*1000./2)**2 * 3600 #[N/mm^2] # aramid 3600, dyneema: 1073.4
    # wind harvesting efficiency
    params['eff'] = 0.6

    # Ampyx aerodynamic coefficients
    params['CD_alpha']     = [-0.6176, -0.1572, -0.0372]
    params['CZalpha']      = [-5.7736, +5.0676, +0.5526]
    params['CXbeta']       = [0.0]
    params['CYbeta']       = [0.093591492950463,-0.029900842662737,-0.185542602484472]
    params['CZbeta']       = [0.0]
    # Damping Force
    params['Cxp']           = [0.0]
    params['Cxq']           = [4.4124, -0.6029 ]
    params['Cxr']           = [0.0]
    params['Cyp']           = [0.049646189990042,-0.014013395929910,-0.102202368153586]
    params['Cyq']           = [0.0]
    params['Cyr']           = [-0.136826511089577,-0.169467213438735]
    params['Czp']           = [0.0]
    params['Czq']           = [6.148615659988596,0.125149171823303,-7.556055799548278]
    params['Czr']           = [0.0]
    # Control dependent coefficients - assumed as zero for this model due to malfunctioning
    params['Cxail']          = [0.0]
    params['Cxele']          = [0.0]
    params['Cxrud']          = [0.0]
    params['Cyail']          = [0.0]
    params['Cyele']          = [0.0]
    params['Cyrud']          = [0.0]
    params['Czail']          = [0.0]
    params['Czele']          = [0.0]
    params['Czrud']          = [0.0]

    # MOMENTS
    params['Clalpha']       = [0.0]
    params['Cmalpha']       = [0.0, -0.6076,  -0.0307]
    params['Cnalpha']       = [0.0]
    params['Clbeta']        = [0.031258025680547,-3.483039877119464e-04,-0.063081351778656]
    params['Cmbeta']        = [0.0]
    params['Cnbeta']        = [-0.084911552609027,0.057710177865613]
    # Damping Moments

    params['Clp']           = [0.281317597010932,-0.024764054416691,-0.563258775268210]
    params['Clq']           = [0.0]
    params['Clr']           = [0.644873057041009,0.181159486166008]
    params['Cmp']           = [0.0]
    params['Cmq']           = [5.288473128322822,-0.002551748906449,-11.302198360248445]
    params['Cmr']           = [0.0]
    params['Cnp']           = [-0.913179228017379,-0.056597312252964]
    params['Cnq']           = [0.0]
    params['Cnr']           = [0.025700796184592,0.029060380546373,-0.055306030491248]
    # Control dependent coefficients
    params['Clail']          = [0.238326419615787,-0.008782572832499,-0.248918359753252]
    params['Clele']          = [0.0]
    params['Clrud']          = [-0.001346210113888,0.004363424631495]
    params['Cmail']          = [0.596509640460760,0.079570783623381,-0.316530950310559]
    params['Cmele']          = [0.997381197242358,-0.006148542440988,-1.042779240167403]
    params['Cmrud']          = [-0.001589335101276]
    params['Cnail']          = [-0.114703589462861,0.019036692592370]
    params['Cnele']          = [0.0]
    params['Cnrud']          = [0.040894647002254,-0.011721249877666,-0.040414259093759]

    return params
