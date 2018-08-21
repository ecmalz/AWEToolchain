'''
MAIN Code to run the optimisation of a drag kite system
based on direct Collocation.py
winddata: polynomials direct from wind data are used to represent the wind field.

--> CHOOSE if FULL homotopy should be solved or starting by running through WIND profiles

Python Version 2.7 / Casadi version 3.3.0
- Author: Elena Malz, Chalmers 2017
'''

kitesize = '2'
import sys
sys.path.append(r"/usr/local/casadi-py27-v3.3.0/")
sys.path.append('../')
import casadi as ca
import casadi.tools as ca
import numpy as np
import pickle
import _OCPmodel as OCP
import _fcts
from _aero_poly import aero_drag_ampyx as aero
sys.path.append('/Users/elenama/Documents/Git/Dragonfly3.0')

if kitesize == '0.666ov':    from parameters import drag_666kW_oversized as initial_params; print 'kitesize 0.666MW oversized is chosen'
elif kitesize == '0.666':   from parameters import drag_666kW_Makani as initial_params; print 'kitesize 0.666MW is chosen'
elif kitesize == '2':       from parameters import drag_2MW as initial_params; print 'kitesize 2MW is chosen'
else:  print 'chosen kitesize not available'

from Collocation import collocate
# ---  Import wind data (all have to be same size) ( from plot Clusters)
with open('../../../Data/polycoeff.dat', 'r') as f:
    (taux_opt,tauy_opt,altitude) = pickle.load(f)

k_coeff = taux_opt.shape[0]


# -----------------------------------------------------------
params     = initial_params()
nk         = 20     # Control discretization     % for a longer time horizon nk has to be high
d          = 3      # Degree of interpolating polynomial yields integration order 5
tf_init    = 9.0    # End time initial guess

# ------------------------
# MODEL
# ------------------------
print '#################################'
print '###        Building OCP...    ###'
print '#################################'

Kite = OCP.Drag(nk,d,k_coeff) # Get model parameters
xd = Kite.xd
u   = Kite.u
xa  = Kite.xa
p   = Kite.p
l   = Kite.l
puni = Kite.puni

out, out_fun = Kite.createoutputs(aero, params)
# -------------
# AERODYNAMICS
# -------------
dynamics = Kite.dynamics(aero,params)  # Create implicit function of system dynamics

# -----------------------------
# DISCRETIZATION / SET UP NLP
# -----------------------------
V, P, coll_cstr, continuity_cstr, Output = collocate(xd,xa,u,p,puni,l,nk,d,dynamics, out_fun,out)
Out_fun = ca.Function('Ofcn',[V,P],[Output])
def get_Output(V,P):
    return Output(Out_fun(V,P))

# --------------------------------------------------
# ADD PATH CONSTRAINTS AND BOUNDARY CONDITIONS
# --------------------------------------------------
tether_cstr, DCM_cstr, periodic_cstr = Kite.PathConstraints(V,params)

# --- OUTPUT CONSTRAINTS ----
output4cstr = get_Output(V,P)

E_cstr          = []  # Initial Energy to zero
AoA_cstr        = [] # limit angle of attack
sslip_cstr      = [] # limit side slip
F_tether_cstr   = [] # maximal allowed tether force
RatedP_cstr     = []

AoA_cstr.append(output4cstr['AoA_deg'])
sslip_cstr.append(output4cstr['sslip_deg'])
E_cstr.append( V['Xd',0,0, 'E'] )
F_tether_cstr.append(output4cstr['tethercstr'])
RatedP_cstr.append(V['Xd',-1,-1,'E']/float(params['ScalePower'])/V['tf'])

# --- Struct of all constraints:
g = ca.struct_MX(
              [
               ca.entry('collocation',     expr=coll_cstr),
               ca.entry('continuity',      expr=continuity_cstr),
               ca.entry('DCM orthogonality', expr=DCM_cstr),
               ca.entry('periodicity',     expr=periodic_cstr),
               ca.entry('tether',          expr = tether_cstr),
               ca.entry('E_cstr',          expr = E_cstr)
            ]
              )

h = ca.struct_MX([
               ca.entry('AoA_cstr',        expr = AoA_cstr),
               ca.entry('sslip_cstr',      expr = sslip_cstr),
               ca.entry('F_tether_cstr',   expr = F_tether_cstr),
               ca.entry('max_power',      expr = RatedP_cstr)
            ])

cstr_struct = ca.struct_MX([ca.entry('g', expr = g, struct = g),
                         ca.entry('h', expr = h,  struct = h)])

gh = cstr_struct(ca.vertcat(g,h))

# --- Concatenate constraints
lbg = cstr_struct()
ubg = cstr_struct()
lbg['h','AoA_cstr'] = -15  # in degrees
ubg['h','AoA_cstr'] =  15
lbg['h','sslip_cstr'] = -30
ubg['h','sslip_cstr'] =  30
lbg['h','F_tether_cstr'] = -ca.inf
ubg['h','F_tether_cstr'] = 0.
lbg['h','max_power'] = -params['Pnom'] * 1e6 #  [W]
ubg['h','max_power'] = ca.inf

# --------------------------------
# OBJECTOVE FUNCTION
# ------------------------------
Cost = Kite.CostFunction(out,V,P, params)


print '##################################'
print '###         Initializing...    ###'
print '##################################'
# --------------
# INITIALIZE STATES
# -------------
vars_init = V()
tau_roots = ca.collocation_points(d, 'radau')
tau_roots = ca.veccat(0, tau_roots)

for k in range(nk):
    for j in range(d+1):
        t = (k + tau_roots[j])*tf_init/float(nk)
        guess = Kite.initial_guess(t,tf_init,params)

        vars_init['Xd',k,j,'q']  = guess['q']
        vars_init['Xd',k,j,'dq'] = guess['dq']
        vars_init['Xd',k,j,'w']  = guess['w']
        vars_init['Xd',k,j,'R']  = guess['R']

# --------------
# BOUNDS
# -------------
vars_lb   = V(-ca.inf)
vars_ub   = V(ca.inf)

vars_lb['Xd',:,:,'q',2] = 1.5 * params['bref']
vars_lb['Xd',:,:,'coeff'] = -25*ca.pi/180. #-------------------------------NOTE THIS IS CHANGED!!!!
vars_ub['Xd',:,:,'coeff'] =  25*ca.pi/180.
# vars_ub['Xd',:,:,'q',-1] =  800. ############# IF NO TETHER DRAG
vars_lb['XA',:]           = 0   # tension should be positive
vars_init['l']            = params['l']
vars_lb['l']              = params['l']
vars_ub['l']              = params['l']


# ---------------------
# PARAMETER VALUES
# --------------------

# --- Build reference parameters, references of cost function should match the initial guess
p_num = P()

for name in xd.keys():
    p_num['p',:,:,'ref',name] = vars_init['Xd',:,:,name]

# --- wind velocities
p_num['puni','weights','q']         = 1.
p_num['puni','weights','dq']        = 1.
p_num['puni','weights','R']         = 1.
p_num['puni','weights','w']         = 1.
p_num['puni','weights','coeff']     = 0.01

p_num['puni','weights','u']         = 0.01
p_num['puni','weights','T']         = 0.01
p_num['puni','weights','dDrag']     = 0.001
p_num['puni','weights','dcoeff']    = 10.
p_num['puni','weights','AoA']       = 1.
p_num['puni','weights','sslip']     = 1.

# --- Wind data (polynomial or log function)
p_num['puni','tau_x']               = taux_opt
p_num['puni','tau_y']               = tauy_opt
p_num['puni','altitude']            = altitude
# p_num['puni','wind']                = params['wind0']
# p_num['puni',:,:,'z0']   = 0.15
p_num['tf'] = tf_init
vars_init['tf'] =  vars_lb['tf'] = vars_ub['tf']  = p_num['tf']


## --------------------
## SOLVE THE NLP
## --------------------

# --- Allocate an NLP solver
nlp = {'x': V, 'p': P, 'f': Cost, 'g': gh}
# Set options
opts = {}
opts["expand"] = True
opts['ipopt.linear_solver'] = 'ma57'
opts['ipopt.max_iter']  = 1000
opts['ipopt.tol']       = 1e-8

solver = ca.nlpsol("solver", "ipopt", nlp, opts)






#########################

def create_initarg(w,free_l,vars_init):
    p_num['puni','wind']         = w

    arg = {}
    # Bounds on x
    arg['lbx'] = vars_lb
    arg['ubx'] = vars_ub
    # Bounds on g
    arg['lbg'] = lbg
    arg['ubg'] = ubg
    Homotopy_step = 0.1

    for gamma_value in list(np.arange(0,1.+Homotopy_step,Homotopy_step)):

        p_num['puni','gam'] = gamma_value
        p_num['toggle_to_energy'] = 0.
        external_extra_text = ['HOMOTOPY SOLVE FOR  GAMMA %.1f' % gamma_value]
        # Initial condition
        arg['x0']   = vars_init
        arg['p']   = p_num   # hand over the parameters to the solver

        # Solve the problem
        print '   '
        print 'Solve for gamma:   ',p_num['puni','gam'] #PARAMETER value for homotopy
        print '   '
        res = solver(**arg)
        stats = solver.stats()
        # print stats['t_mainloop.proc']
        # print stats['t_mainloop.wall']
        assert stats['return_status'] in ['Solve_Succeeded']
        print '   '
        print 'Solved for gamma:  ',p_num['puni','gam'] #PARAMETER value for homotopy
        print '   '

        # -----------------------------------------------------------------------------
        # UPDATE INITIAL GUESS
        # -----------------------------------------------------------------------------
        arg['lam_x0'] = res['lam_x']
        vars_init = V(res['x'])

    init_opt = vars_init
    # -----------------------------------------------------------------------------
    # POWER OPTIMISATION
    # -----------------------------------------------------------------------------
    # Using homopty for changing cost function.
    # Shifting from tracking to power optimisation
    print "#####################################################"
    print "#####################################################"
    print "#####################################################"
    print "#########                                   #########"
    print "#########    STARTING POWER OPTIMIZATION    #########"
    print "#########                                   #########"
    print "#####################################################"
    print "#####################################################"
    print "#####################################################"

    # --- Ensure no fictitious forces
    p_num['puni','gam'] = 1.

    #tf_iterations = [] #record orbit periods
    Homotopy_step = 0.1
    toggle_table = list(np.arange(Homotopy_step,1.+ Homotopy_step,Homotopy_step))

    for toggle_value in toggle_table:

        # --- Attribute previous guess
        arg['x0'] = vars_init

        # --- Update toggle value
        p_num['toggle_to_energy'] = toggle_value
        arg['p']  = p_num
        external_extra_text = ['POWER OPTIMIZATION; TOGGLE %.1f - FIXED TIME' % toggle_value]
        # --- Solve the problem
        print "Solve for toggle =", toggle_value
        res = solver(**arg)

        # --- Retrieve the solution, re-assign as new guess
        arg['x0']            = res['x']
        arg['lam_x0']        = res['lam_x']
        # p_num['tf_previous'] = V(res['x'])['tf']
        arg['p']             = p_num

        # --- Report some stuff...
        print "Solved for toggle =", toggle_value, " Period = ", float(V(res['x'])['tf'])
        vars_init = V(res['x'])


    opt = V(res['x'])

    print "\n\n\n"
    print "Average Power = ", -opt['Xd',-1,-1,'E']/float(params['ScalePower'])/opt['tf'], "  Orbit period = ", opt['tf']

    if free_l:
        print "#####################################################"
        print "#####################################################"
        print "#####################################################"
        print "#########                                ############"
        print '#########     OPTIMISE TETHER LENGTH     ###########'
        print "#########                                ############"
        print "#####################################################"
        print "#####################################################"
        print "#####################################################"

        vars_init         = opt
        vars_lb['l']      = 50.
        vars_ub['l']      = 10000.
        arg['x0']         = vars_init
        arg['lbx']        = vars_lb
        arg['ubx']        = vars_ub
        res               = solver(**arg)
        print "\n\n\n"
        print "Average Power = ", -opt['Xd',-1,-1,'E']/float(params['ScalePower'])/opt['tf'], "  Orbit period = ", opt['tf']
        print 'optimal length', V(res['x'])['l']

    print "#####################################################"
    print "#####################################################"
    print "#####################################################"
    print "#########                                   #########"
    print "#########          OPEN FINAL TIME          #########"
    print "#########                                   #########"
    print "#####################################################"
    print "#####################################################"
    print "#####################################################"

    vars_lb['tf'] = 0.
    vars_ub['tf'] = ca.inf

    arg['lbx'] = vars_lb
    arg['ubx'] = vars_ub
    arg['x0'] = res['x']


    external_extra_text = ['RELEASE TIME - final solve']
    res = solver(**arg)

    #------------------------------
    # RECEIVE SOLUTION  & SAVE DATA
    #------------------------------
    opt = V(res['x'])

    outputs = Output(Out_fun(V,P))
    val_init = get_Output(vars_init, p_num)
    val_opt = get_Output(opt,p_num)

    # with open('../../../Data/drag_dat/solution.dat','w') as f:
    #     pickle.dump((val_opt, opt, nk, d),f)

    # with open('../../../Data/drag_dat/init.dat','w') as f:
    #     pickle.dump((val_init, vars_init, nk, d),f)
    #

    with open('arg_Pcap' +str(p_num['puni','wind'])+ '.dat','w') as f:
        pickle.dump((arg, res, p_num),f)


    print '##############################'
    print 'saved:   arg_Pcap' +str(p_num['puni','wind'])+ '.dat'
    print '##############################'
    #print val_opt['speed']
    #print val_opt['windspeed_shear']
    # E_final   = Energy_Cost_fun(opt,p_num)
    # Lifting   = lift_Cost_fun(opt)
    # Tracking  = Tracking_Cost_fun(opt,p_num)
    # Cost      = totCost_fun(opt,p_num)
    # Reg       = Reg_Cost_fun(opt,p_num)
    #
    # with open('../../../Data/drag_dat/cost.dat', 'w') as f:
    #     pickle.dump((E_final, Lifting, Tracking, Cost, Reg), f)

    # --------------------------------------
    # PRINT OUT ....
    # --------------------------------------
    print "\n\n\n\n\n\n"
    print 'wind', p_num['puni','wind']
    print 'tf', tf_init
    print "Average Power = ", -opt['Xd',-1,-1,'E']/float(params['ScalePower'])/opt['tf'], "[W]  Orbit period = ", opt['tf']
    return opt, val_opt
