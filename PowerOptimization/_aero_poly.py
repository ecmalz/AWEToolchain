
'''
simple AERODYNAMIC MODEL of AWE system (pumping mode)
### creates _POLYNOMIAL_ for wind profile based on different wind data
takes states and inputs and creates aerodynamic forces and moments
dependent on the position of the kite.
Aerodynamic coefficients are assumptions.
Python Version 2.7 / Casadi version 3.3.0
- Author: Elena Malz, Chalmers 2017
'''

import casadi as ca
import casadi.tools as ca
import sys


import numpy as np
from numpy import pi
sys.path.append('..')
import fcts


def aero(xd, xa, p, puni,l, params):
    [q, dq, R, w, coeff, E, Drag ] = xd[...]
    # --- Build current wind shear profile
    #
    Lagrange_polynomial_x = fcts.Lagrange_poly(puni['altitude'],puni['tau_x'])
    Lagrange_polynomial_y = fcts.Lagrange_poly(puni['altitude'],puni['tau_y'])
    # --- wind vector & apparent velocity
    xwind = Lagrange_polynomial_x(q[-1])
    ywind = Lagrange_polynomial_y(q[-1])
    v_app = dq - np.array([xwind,ywind,0])


    # xwind = puni['wind']*(q[-1]/params['windShearRefAltitude'])**0.15   #NOTE:AUTOMATE THAT!!!!!!
    # v_app = dq - np.array([xwind,0,0])


    # --- calculate angle of attack and side slip (convert apparent velocity to body frame)
    AoA       = -ca.mtimes(R[6:9].T,(v_app))/ca.mtimes(R[0:3].T,(v_app))
    sslip     = ca.mtimes(R[3:6].T,(v_app))/ca.mtimes(R[0:3].T,(v_app))

    params['rho'] = fcts.rho(q[-1])
    CF, CM, Cl, Cm, Cn = aero_coeffs(AoA, sslip, v_app, w,  coeff , params)

    F_aero = 0.5 * params['rho'] * params['sref'] * ca.norm_2(v_app) * (ca.cross(CF[2]*v_app,R[3:6]) + CF[0]*(v_app) )
    reference_lengths = ca.vertcat([params['bref'], params['cref'], params['bref']])
    M_aero = 0.5 * params['rho'] * params['sref'] *  ca.norm_2(v_app)**2 * CM  * reference_lengths

    F_drag      = -Drag * R[0:3]*params['eff']   # Drag in opposite x-direction of kite (Props)
    F_tether    = q*xa             # Force of the tether at the kite
    # m           = params['mK'] + 1./3*params['mT']
    m   = params['mK'] + 1./3*params['tether_density']*l    #mass of kite and tether

    F_gravity   = m*params['g']

    # --- Tether Drag, defined in + x direction , therefore minus sign
    C_tet        = 0.4
    Tether_drag  =  -1./8 * params['rho'] * params['tether_diameter'] * C_tet * l * ca.norm_2(v_app)**2 * R[0:3]
    # Tether_drag  =  -  1./6 * params['rho'] * params['tether_diameter'] * C_tet * l * ca.norm_2(v_app)**2 * R[0:3]

    outputs = {}
    outputs['v_app']            = v_app
    # outputs['rho']              = params['rho']
    outputs['speed']            = ca.norm_2(v_app)
    outputs['windspeed_shear']  = xwind
    # outputs['windspeedy_shear'] = ywind
    outputs['momentconst']      = ([Cl,Cm,Cn])
    outputs['AoA']              = AoA
    outputs['sslip']            = sslip
    outputs['AoA_deg']          = AoA  * 180/pi
    outputs['sslip_deg']        = sslip * 180/pi
    outputs['CL']               = CF[2]
    outputs['CD']               = -CF[0]
    outputs['F_aero']           = F_aero
    outputs['F_drag']           = F_drag
    outputs['F_tether_scaled']  = F_tether   # This is scaled so Force/kg
    outputs['F_tether']         = F_tether * m
    outputs['tethercstr']       = (xa * l * m) - pi/3. * (params['tether_diameter']*1000/2.)**2 * 3600   # tdiameter in m, l*xa = F_tether instead of norm (q*xa)
    # Tensile strength of aramid: 3600N/mm^2, while dyneema has only 1700?
    outputs['M_aero']           = M_aero
    outputs['power']            = ca.mtimes(F_drag.T,dq)
    outputs['Tether_drag']      = Tether_drag
    return (F_aero, M_aero, F_tether, F_drag, F_gravity, Tether_drag), outputs

def aero_coeffs(alpha, beta, v_app, omega,  phi , params):
    ''' gives all aerodynamic coefficients '''

    alphaDeg = alpha * 180./np.pi
    CL = params['CL0'] + params['CLalpha'] * alphaDeg
    CD = params['CD0'] + params['CDalpha'] * alphaDeg + params['CDalphaSq'] * alphaDeg**2.0 + 2.*beta**2 + 0.2*alpha**2
    # CL = (params['CL0']-0.05) + (params['CLalpha']+0.03) * alphaDeg
    # CD = params['CD0']+0.002 + (params['CDalpha']+0.002) * alphaDeg + (params['CDalphaSq']-0.0003) * alphaDeg**2.0 + 2.*beta**2 + 0.2*alpha**2


    CFx_0 = - CD   # theoretical adding rotor drag but no rotors
    CFy_0 = 0.     #2*pi*0.1 * beta            #num_pylons * CY_b * pylon_sref/params['sref']
    CFz_0 = CL

    CMx_0 = 0.
    CMy_0 = 0.
    CMz_0 = 0.

    # pqr - DAMPING

    CFx_pqr = 0.
    CFy_pqr = 0.
    CFz_pqr = 0.

    omega_hat = omega/(2.*ca.norm_2(v_app))
    omega_hat[0] *= params['bref']
    omega_hat[1] *= params['cref']
    omega_hat[2] *= params['bref']

    p = omega_hat[0]
    q = omega_hat[1]
    r = omega_hat[2]

    # roll moment, about ehat1
    Cl = params['Clp'] * p + params['Clr'] * r + params['Clbeta'] * beta
    # pitch moment, about ehat2
    Cm = params['Cmq'] * q + params['Cmalpha'] * alpha
    # yaw moment, about ehat3
    Cn = params['Cnp'] * p + params['Cnr'] * r  + params['Cnbeta'] * beta

    CMx_pqr = Cl
    CMy_pqr = Cm
    CMz_pqr = Cn

    # surfaces
    CFx_surfs = 0.
    CFy_surfs = 0.
    CFz_surfs = 0.
    CMx_surfs = phi[0] #* 0.1
    CMy_surfs = phi[1] #* 0.1
    CMz_surfs = phi[2] #* 0.1


    CFx = CFx_0 + CFx_pqr + CFx_surfs
    CFy = CFy_0 + CFy_pqr + CFy_surfs
    CFz = CFz_0 + CFz_pqr + CFz_surfs
    CMx = CMx_0 + CMx_pqr + CMx_surfs
    CMy = CMy_0 + CMy_pqr + CMy_surfs
    CMz = CMz_0 + CMz_pqr + CMz_surfs

    CF_wind = ca.vertcat(CFx, CFy, CFz) # in fixed frame
    CM_cad = ca.vertcat(CMx, CMy, CMz)  # in body frame

    return CF_wind, CM_cad, Cl, Cm, Cn

def aero_drag_ampyx(xd, xa, p, puni, l, params):
    [q, dq, R, w, coeff, E, Drag ] = xd[...]
    # --- Build current wind shear profile
    Lagrange_polynomial_x = fcts.Lagrange_poly(puni['altitude'],puni['tau_x'])
    Lagrange_polynomial_y = fcts.Lagrange_poly(puni['altitude'],puni['tau_y'])
    # --- wind vector & apparent velocity
    xwind = Lagrange_polynomial_x(q[-1])
    ywind = Lagrange_polynomial_y(q[-1])
    v_app = dq - np.array([xwind,ywind,0])


    # xwind = puni['wind']*(q[-1]/params['windShearRefAltitude'])**0.15   #NOTE:AUTOMATE THAT!!!!!!
    # v_app = dq - np.array([xwind,0,0])


    v_app_body      = ca.mtimes(R.T,v_app)

    # calculate angle of attack and side slip (convert apparent velocity to body frame)
    AoA       = -ca.mtimes(R[6:9].T,(v_app))/ca.mtimes(R[0:3].T,(v_app))
    sslip     = ca.mtimes(R[3:6].T,(v_app))/ca.mtimes(R[0:3].T,(v_app))

    CF, CM, CD, CL = aero_coeffs_ampyx(AoA, sslip, v_app, w,  coeff , params)

    W0 = v_app / ca.norm_2(v_app)
    W2 = ca.cross(v_app, R[3:6]) / ca.norm_2(v_app)
    W1 = ca.cross(W2,W0)
    W = ca.vertcat(W0.T, W1.T, W2.T).T

    #NED to zUp
    #NED2zup =  [ [1, 0, 0], [0, -1, 0], [0, 0, -1]]
    #NED2zup = np.reshape(NED2zup, (3,3))
    #CF      = ca.mtimes(NED2zup,CF)
    F_aero  = ca.mtimes(W,CF)
    M_aero  = CM
    # F_aero = 0.5 * params['rho'] * params['sref'] * ca.norm_2(v_app) * (ca.cross(-CF[2]*v_app,R[3:6]) + CF[0]*(v_app) )
    # reference_lengths = ca.vertcat([params['bref'], params['cref'], params['bref']])
    # M_aero = 0.5 * params['rho'] * params['sref'] *  ca.norm_2(v_app)**2 * CM # * reference_lengths

    F_drag   = -Drag * R[0:3]* params['eff']   # Drag in opposite x-direction of kite (Props)
    F_tether = q*xa             # Force of the tether at the kite
    m = params['mK'] +  1./3*params['tether_density']*l
    # m   = params['mK'] + 1./3*params['tether_density']*ltet    #mass of kite and tether

    F_gravity = m*params['g']

    # Tether Drag
    # defined in + x direction , therefore minus sign
    C_tet        = 1.2
    # Tether_drag  =  - 1./8 * params['rho'] * params['tether_diameter'] * C_tet * ltet * ca.norm_2(v_app)**2 * ca.veccat(ca.cos(AoA), 0, ca.sin(AoA))
    Tether_drag  =  - 1./8 * params['rho'] * params['tether_diameter'] * C_tet * l * ca.norm_2(v_app)* v_app

    Tether_drag_body = Tether_drag               #!!!!----NOTE: It should be nav frame from the beginning, since we use v_app_nav, right?

    outputs = {}
    outputs['v_app']    = v_app
    outputs['speed']    = ca.norm_2(v_app)
    outputs['windspeed_shear'] = xwind
    outputs['AoA']      = AoA
    outputs['sslip']    = sslip
    outputs['AoA_deg']  = AoA  * 180/pi
    outputs['sslip_deg']= sslip * 180/pi
    outputs['CL']       = CL[0][0]
    outputs['CD']       = -CD[0][0]
    outputs['F_aero']   = F_aero
    outputs['F_aero_b'] = ca.mtimes(R.T,F_aero)
    outputs['F_drag']   = F_drag
    outputs['F_tether_scaled'] = F_tether   # This is scaled so Force/kg
    outputs['F_tether'] = F_tether*m
    outputs['tethercstr'] = (xa * l * m) - pi/3. * (params['tether_diameter']*1000/2.)**2 * 3600   # tdiameter in m, l*xa = F_tether instead of norm (q*xa)
    outputs['M_aero']   = M_aero
    outputs['power']    = ca.mtimes(F_drag.T,dq)
    outputs['Tether_drag'] = Tether_drag
    return (F_aero, M_aero, F_tether, F_drag, F_gravity, Tether_drag), outputs

def aero_coeffs_ampyx(alpha, beta, v_app, omega,  phi , params):
    """ loading coefficients from ampyx AP2 and calculateing F and M with Ampyx style (NED frame) """

    CZ_alpha = ca.polyval(params['CZalpha'],alpha)
    CX_alpha = ca.polyval(params['CD_alpha'],alpha)

    CX_beta = ca.polyval(params['CXbeta'],alpha)
    CX_p    = ca.polyval(params['Cxp'],alpha)
    CX_q    = ca.polyval(params['Cxq'],alpha)
    CX_r    = ca.polyval(params['Cxr'],alpha)
    CX_deltaA  = ca.polyval(params['Cxail'],alpha)
    CX_deltaE  = ca.polyval(params['Cxele'],alpha)
    CX_deltaR  = ca.polyval(params['Cxrud'],alpha)

    CY_beta = ca.polyval(params['CYbeta'],alpha)
    CY_p    = ca.polyval(params['Cyp'],alpha)
    CY_q    = ca.polyval(params['Cyq'],alpha)
    CY_r    = ca.polyval(params['Cyr'],alpha)
    CY_deltaA  = ca.polyval(params['Cyail'],alpha)
    CY_deltaE  = ca.polyval(params['Cyele'],alpha)
    CY_deltaR  = ca.polyval(params['Cyrud'],alpha)

    CZ_beta = ca.polyval(params['CZbeta'],alpha)
    CZ_p    = ca.polyval(params['Czp'],alpha)
    CZ_q    = ca.polyval(params['Czq'],alpha)
    CZ_r    = ca.polyval(params['Czr'],alpha)
    CZ_deltaA  = ca.polyval(params['Czail'],alpha)
    CZ_deltaE  = ca.polyval(params['Czele'],alpha)
    CZ_deltaR  = ca.polyval(params['Czrud'],alpha)

    Cl_beta = ca.polyval(params['Clbeta'],alpha)
    Cl_p    = ca.polyval(params['Clp'],alpha)
    Cl_q    = ca.polyval(params['Clq'],alpha)
    Cl_r    = ca.polyval(params['Clr'],alpha)
    Cl_deltaA  = ca.polyval(params['Clail'],alpha)
    Cl_deltaE  = ca.polyval(params['Clele'],alpha)
    Cl_deltaR  = ca.polyval(params['Clrud'],alpha)

    Cm_alpha = ca.polyval(params['Cmalpha'],alpha)
    Cm_beta = ca.polyval(params['Cmbeta'],alpha)
    Cm_p    = ca.polyval(params['Cmp'],alpha)
    Cm_q    = ca.polyval(params['Cmq'],alpha)
    Cm_r    = ca.polyval(params['Cmr'],alpha)
    Cm_deltaA  = ca.polyval(params['Cmail'],alpha)
    Cm_deltaE  = ca.polyval(params['Cmele'],alpha)
    Cm_deltaR  = ca.polyval(params['Cmrud'],alpha)

    Cn_beta = ca.polyval(params['Cnbeta'],alpha)
    Cn_p    = ca.polyval(params['Cnp'],alpha)
    Cn_q    = ca.polyval(params['Cnq'],alpha)
    Cn_r    = ca.polyval(params['Cnr'],alpha)
    Cn_deltaA  = ca.polyval(params['Cnail'],alpha)
    Cn_deltaE  = ca.polyval(params['Cnele'],alpha)
    Cn_deltaR  = ca.polyval(params['Cnrud'],alpha)



    # --- pqr - DAMPING ---
    phat = omega[0]/(2.*ca.norm_2(v_app)) * params['bref']
    qhat = omega[1]/(2.*ca.norm_2(v_app)) * params['cref']
    rhat = omega[2]/(2.*ca.norm_2(v_app)) * params['bref']

    CXw   = CX_beta * beta + CX_p * phat + CX_q * qhat + CX_r * rhat + CX_alpha
    CXphi = CX_deltaE * phi[1] + CX_deltaA * phi[0] + CX_deltaR * phi[2]
    CYw   = CY_beta * beta + CY_p * phat + CY_q * qhat + CY_r * rhat
    CYphi = CY_deltaE * phi[1] + CY_deltaA * phi[0] + CY_deltaR * phi[2]
    CZw   = CZ_beta * beta + CZ_p * phat + CZ_q * qhat + CZ_r * rhat + CZ_alpha
    CZphi = CZ_deltaE * phi[1] + CZ_deltaA * phi[0] + CZ_deltaR * phi[2]

    CX = CXw + CXphi
    CY = CYw + CYphi
    CZ = CZw + CZphi

    Cl = Cl_beta * beta + Cl_p * phat + Cl_q * qhat + Cl_r * rhat + \
        Cl_deltaE * phi[1] + Cl_deltaA * phi[0] + Cl_deltaR * phi[2]
    Cm = Cm_beta * beta + Cm_p * phat + Cm_q * qhat + Cm_r * rhat + Cm_alpha +\
        Cm_deltaE * phi[1] + Cm_deltaA * phi[0] + Cm_deltaR * phi[2]
    Cn = Cn_beta * beta + Cn_p * phat + Cn_q * qhat + Cn_r * rhat + \
        Cn_deltaE * phi[1] + Cn_deltaA * phi[0] + Cn_deltaR * phi[2]

    qbar = ca.norm_2(v_app)**2 * 0.5 * params['rho']

    FX                     = CX * qbar * params['sref']
    FY                     = CY * qbar * params['sref']
    FZ                     = CZ * qbar * params['sref']
    F_AeroBody             = ca.vertcat(FX, FY, FZ)

    MX                     = qbar * params['sref'] * params['bref'] * Cl
    MY                     = qbar * params['sref'] * params['cref'] * Cm
    MZ                     = qbar * params['sref'] * params['bref'] * Cn
    M_Aero                 =  ca.vertcat(MX, MY, MZ)
    return F_AeroBody,M_Aero, CX_alpha, CZ_alpha
