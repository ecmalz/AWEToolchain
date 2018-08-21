import sys
sys.path.append(r"/usr/local/casadi-py27-v3.3.0/")
import casadi as ca
import casadi.tools as ca
# from parameters import AP2
import numpy as np




def skew(a,b,c):
    " creates skew-symmetric matrix"
    d =  ca.blockcat([[0.,        -c, b],
                   [c,  0.,       -a],
                   [-b,  a,    0.]])
    return d
#

class Drag(object):
    '''model for pumping mode with right hand side coordinate, x down'''

    def __init__(self,nk=20,d=3,k_coeff = 10):

        self.nk = nk
        self.d = d

        self.l   = ca.SX.sym('l')
        self.xa  = ca.SX.sym('xa')      # algebraic state
        self.p   = ca.SX.sym('p')       # parameters


        self.xd = ca.struct_symSX([(
                            ca.entry('q', shape = (3,1)),      # position
                            ca.entry('dq', shape = (3,1)),     # velocity of the kite
                            ca.entry('R', shape = (3,3)),      # rotation matrix DCM: transfors from body to fixed frame
                            ca.entry('w', shape = (3,1)),      # angular velocity
                            ca.entry('coeff', shape = (3,1)),  # moment coefficients Cr, Cp, Cy
                            ca.entry('E'),                     # generated energy
                            ca.entry('Drag'),                  # Drag Force
                            )])


        self.u = ca.struct_symSX([(                                 # 3 dimensional control input + torque
                           ca.entry('u', shape = (3,1)),
                           ca.entry('T', shape = (3,1)),
                           ca.entry('dcoeff', shape = (3,1)),
                           ca.entry('dDrag'),                     # control input (change in drag)
                           )])



        ref = ca.struct_symSX( [ ca.entry(name, shape = self.xd[name].shape) for name in self.xd.keys() ] + \
                            [ ca.entry(name, shape =  self.u[name].shape) for name in  self.u.keys() ]
                             )

        self.weights =  ca.struct_symSX(
                                [ca.entry(name, shape = self.xd[name].shape) for name in self.xd.keys()]  + \
                                [ca.entry(name, shape = self.u[name].shape) for name in  self.u.keys()]  + \
                                [ca.entry('AoA'),
                                ca.entry('sslip'),
                                ca.entry('l')
                                ])

        self.p = ca.struct_symSX([(
                           ca.entry('ref', struct = ref),                             # tracking reference for cost function
                                )])


        self.puni = ca.struct_symSX([(
                   ca.entry('tau_x',          shape =k_coeff ),               # x parameters for the wind polynomial'
                   ca.entry('tau_y',          shape =k_coeff ),               # y parameters for the wind polynomial'
                   ca.entry('altitude',       shape =k_coeff ),               # heights (xdata) for the wind polynomial
                   ca.entry('weights',        struct = self.weights),              # weights for cost function
                   ca.entry('wind'),                                          # base wind (if not polynomial)
                   ca.entry('gam'),                                           # homopty variable
                   )])


    def createoutputs(self,aero,params):
        _,outputs = aero(self.xd, self.xa, self.p, self.puni, self.l, params)
        [q, dq, R, w, coeff, E, Drag ] = self.xd[...]
        outputs['c']  = ( ca.sum1(q**2) - self.l**2 )
        outputs['dc'] = ca.sum1(q*dq)
        out = ca.struct_SX( [ca.entry(name, expr=outputs[name]) for name in outputs] )
        out_fun = ca.Function('outputs', [self.xd , self.xa, self.p, self.puni, self.l], [out])
        return out, out_fun


    def dynamics(self, aero, params):
        ScalePower = params['ScalePower']
        scale = 1000.
        xd  = self.xd
        u   = self.u
        xa  = self.xa
        p   = self.p
        puni = self.puni
        l = self.l
        m   = params['mK'] + 1./3*params['mT']    #mass of kite and tether
        g   = params['g']               # gravity
        A   = params['sref']            # area of Kite
        J   = params['J']               # Kite Inertia

        xddot = ca.struct_symSX([ ca.entry('d'+name, shape = xd[name].shape) for name in xd.keys() ] )
        [q, dq, R, w, coeff, E, Drag ] = xd[...]

        (Fa, M,  F_tether_scaled, F_drag, F_gravity, Tether_drag), outputs = aero(xd, xa, p, puni, l, params)

        wx = skew(w[0],w[1],w[2]) # creating a skrew matrix of the angular velocities


        # Rdot = R*w (whereas w is in the skrew symmetric form)
        Rconstraint = ca.reshape( xddot['dR'] - ca.mtimes(xd['R'],wx),9,1 ) # J*wdot +w x J*w = T
        TorqueEq = ca.mtimes(J,xddot['dw']) + (ca.cross(w.T,ca.mtimes(J,w).T).T - scale*(1-puni['gam'])*u['T']) - puni['gam']*M
        DragForce = -Drag*R[0:3]

        # ------------------------------------------------------------------------------------
        #  DYNAMICS of the system - Create a structure for the Differential-Algebraic Equation
        # ------------------------------------------------------------------------------------
        res = ca.vertcat(
                       xddot['dq']        - xd['dq'],\
                       xddot['dcoeff']    - u['dcoeff'],\
                       xddot['dDrag']     - u['dDrag'],\
                       m*(xddot['ddq'][0]  + F_tether_scaled[0] - A*(1-puni['gam'])*u['u',0]) -  puni['gam']*Fa[0] - F_drag[0] - Tether_drag[0], \
                       m*(xddot['ddq'][1]  + F_tether_scaled[1] - A*(1-puni['gam'])*u['u',1]) -  puni['gam']*Fa[1] - F_drag[1] - Tether_drag[1], \
                       m*(xddot['ddq'][2]  + F_tether_scaled[2] - A*(1-puni['gam'])*u['u',2]) -  puni['gam']*Fa[2] - F_drag[2] - Tether_drag[2]+ F_gravity   , \
                       xddot['dE'] - ScalePower*ca.mtimes(F_drag.T,dq)
                      )

        res = ca.veccat(res, Rconstraint, TorqueEq)
        res = ca.veccat(res, ca.sum1(xd['q']*xddot['ddq'])+ca.sum1(xd['dq']**2))

        # --- System dynamics function (implicit formulation)
        dynamics = ca.Function('dynamics', [xd,xddot,xa,u,p,puni,l],[res])
        return dynamics





    def initial_guess(self,t,tf_init,params):
        x_guess = self.xd()
        inclination    = 30.0*ca.pi/180.0  # the other angle than the one you're thinking of
        dcmInclination = np.array([[ca.cos(inclination), 0.0, -ca.sin(inclination)],
                                   [                0.0, 1.0,                  0.0],
                                   [ca.sin(inclination), 0.0,  ca.cos(inclination)]])
        dtheta = 2.0 * ca.pi/tf_init
        theta  = t * dtheta
        r      = 0.25 * params['l']

        angle = ca.SX.sym('angle')
        x_cir = ca.sqrt(params['l']**2 - r**2)
        y_cir = r * ca.cos(angle)
        z_cir = r * ca.sin(angle)
        init_pos_fun = ca.Function('init_pos',[angle],[ca.mtimes(dcmInclination, ca.veccat(x_cir, y_cir, z_cir))])
        init_vel_fun = init_pos_fun.jacobian()

        ret = {}
        ret['q']    = init_pos_fun(theta)
        ret['dq']   = init_vel_fun(theta,0) * dtheta
        ret['w']    = ca.veccat(0.0, 0.0, dtheta)

        norm_vel = ca.norm_2(ret['dq'])
        norm_pos = ca.norm_2(ret['q'])

        R0    = ret['dq']/norm_vel
        R2    = ret['q']/norm_pos
        R1    = ca.cross(ret['q']/norm_pos,ret['dq']/norm_vel)
        ret['R'] = ca.vertcat(R0.T, R1.T, R2.T).T
        return ret


    def PathConstraints(self,V,params):
        DCM_cstr        = [] # meet R.T*R - I = 0
        periodic_cstr   = [] # periodic optimisation problem
        tether_cstr     = [] # c = 0;

        # --- ROTATION MATRIX CONSTRAINT
        # rotation matrix DCM has to be orthogonal at every stage. It should be valid R0'R0-I = 0 as well as R0'RN I=0.
        # However the constraints, to get in total only 9 constraints, not 18. (1,2)(1,3)(2,3) to the latter an the rest to the former equation.
        R0R0 = ca.mtimes(V['Xd',0,0,'R'].T,V['Xd',0,0,'R'])   - np.eye(3)
        R0RN = ca.mtimes(V['Xd',0,0,'R'].T,V['Xd',-1,-1,'R']) - np.eye(3)

        for k in [0,3,4,6,7,8]:
            DCM_cstr.append(R0R0[k])

        for k in [1,2,5]:
            DCM_cstr.append(R0RN[k])


        # --- PERIODICITY
        # --- add an constraint so that the inital and the final position is the same.
        # --- hereby initial bounds have to be relaxed so that the optimiser decides by itself where to start.

        xd_names = self.xd.keys()
        for name in set(xd_names)-set(['R','E','q','dq']):
            periodic_cstr.append( V['Xd',0,0,name]-V['Xd',-1,-1,name] )

        periodic_cstr.append( V['Xd',0,0, 'q'] - V['Xd',-1,-1, 'q'] + V['Xd',0,0,'dq']*V['vlift'])
        periodic_cstr.append( V['Xd',0,0,'dq'] - V['Xd',-1,-1,'dq'] + V['Xd',0,0, 'q']*V['vlift'])

        periodic_cstr = ca.veccat(*periodic_cstr)


        # -- TETHER LENGTH --- c = 0 at one time point,  dc = 0 leads to LICQ problems!
        tether_cstr.append(ca.sum1(V['Xd',0,0,'q']**2) - V['l']**2 )
        return tether_cstr, DCM_cstr, periodic_cstr



    def CostFunction(self, out, V, P, params):
        """        # --- tracking dissappears slowly in the cost function and Energy maximising appears. at the final step, cost function
                # --- contains maximising energy, lift, SOSC, and regularisation."""
        u = self.u
        p = self.p
        puni = self.puni
        xd = self.xd
        xa = self.xa
        l = self.l
        Lagrange_Tracking = 0
        Lagrange_Regularisation = 0

        # input regularization
        for name in set(u.keys()):
            Lagrange_Regularisation += puni['weights',name][0]*ca.mtimes((u[name]-p['ref',name]).T,u[name]-p['ref',name])

        Lagrange_Regularisation += puni['weights','AoA']*out['AoA']**2
        Lagrange_Regularisation += puni['weights','sslip']*out['sslip']**2

        # --- Initialization tracking
        for name in set(xd.keys())- set(['R','E','Drag']):
            Lagrange_Tracking += puni['weights',name][0]*ca.mtimes((xd[name]-p['ref',name]).T,xd[name]-p['ref',name])
        for k in range(9):
            Lagrange_Tracking += ca.reshape(puni['weights','R'][0]*ca.mtimes((xd['R']-p['ref','R']).T,xd['R']-p['ref','R']),9,1)[k]


        Lagrange_Tracking       = ca.Function('lagrange_track', [xd,xa,u,p,puni,l],[Lagrange_Tracking])
        Lagrange_Regularisation = ca.Function(  'lagrange_reg', [xd,xa,u,p,puni,l],[Lagrange_Regularisation])


        Tracking       = 0
        Regularisation = 0


        for k in range(self.nk):  # V['XA',k,0] is not same time step as V['Xd',k,0] but same result
            ftrack = Lagrange_Tracking(V['Xd',k,0], V['XA',k,0], V['U',k], P['p',k,0],P['puni'], V['l'])
            Tracking += ftrack

            freg = Lagrange_Regularisation(V['Xd',k,0], V['XA',k,0], V['U',k], P['p',k,0],P['puni'], V['l'])
            Regularisation += freg

        E_final             = 10. * V['Xd',-1,-1,'E']    # for maximising final energy
        Tracking_Cost       = (1-P['toggle_to_energy']) * Tracking          #* 1e-3  # Tracking of initial guess
        Regularisation_Cost = Regularisation                                # Regularisation of inputs
        Lift_Cost           = 0.5*V['vlift']**2 #* 1e2                      # Regularisation of inputs
        Energy_Cost         = P['toggle_to_energy'] * (E_final/params['sref'])/V['tf']
        SOSCFix             = 10. * V['Xd',self.nk/4,0,'q',1]**2

        Cost = 0
        Cost = (Tracking_Cost + Regularisation_Cost + Lift_Cost + SOSCFix)/float(self.nk) + Energy_Cost

        return Cost
