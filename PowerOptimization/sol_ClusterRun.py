'''
COMPARE the different polynomial wind profiles.
Using the SUN class to plot
CHOOSE if you evaluate DRAG or PUMP in the beginning of the code

Python Version 2.7 / Casadi version 3.3.0
- Author: Elena Malz, Chalmers 2017
'''


import sys
sys.path.append(r"/usr/local/casadi-py27-v3.3.0/")
import casadi as ca
import casadi.tools as ca
import numpy as np
import pickle
import matplotlib.pyplot as plt
import pylab as pl
# -- Import other codes
sys.path.append('../')
sys.path.append('../../')
import fcts
plt.rcParams['axes.grid'] = True

plt.close('all')

cities = {'Pamplona': '42.81_-1.65', 'Esbjerg': '55.50_8.50', 'Goteborg': '57.5_11.875', 'Paris': '48.86_2.35',\
    'Tenerife':'28.16_-16.60', 'Galway':'53.27_9.05', 'redsea':'21.73_36.63'}

city = 'Goteborg'
directory = '/Users/elenama/Documents/Git/' + city + '16/'
# prefix = '57.5_11.875_20160101_20161231'
coordinate =  cities[city]
prefix = coordinate + '_20160101_20161231'
power = []

def extract_solution(directory, prefix, Pcap):
    import init_Drag_wind as model   # in order to get the CasADi struct

    # --- Load the solution
    normwind = np.load(directory + prefix+'_features.npy')
    timepoints   = normwind.shape[1]
    sol    = dict(); failed = []
    for t in np.arange(timepoints):
        # --- Collecting the data points, which failed the optimisation
        try:
            with open(directory + Pcap + '/timestep'+ str(t) + '_' + Pcap , 'r') as f:
                (average_power, clu, heights, tau_x, tau_y, res) = pickle.load(f)
        except IOError:
            print ('File:' + directory  + Pcap + '/timestep'+ str(t) + '_' + Pcap + ' is not available.')
            failed.append(t)
            continue

        power.append(float(average_power)*1e-6)
        opt = model.V(res['x'])
        outputs  = model.Output(model.Out_fun(model.V,model.P))
        # val_init = model.get_Output(vars_init,model.p_num)
        val_opt  = model.get_Output(opt,model.p_num)
        sol[t] = [opt,val_opt]

    with open(directory+'power_'+Pcap,'w') as f:
        pickle.dump((power),f)
    with open(directory+Pcap + '/failed_features_'+Pcap,'w') as f:
        pickle.dump((failed),f)
    return sol, power, failed

def load_duration_curve(power):
    sorted(power)
    return  sorted(power)[-1:0:-1]

def operation_altitudes(sol):
    maxheight = [max(sol[k][0]['Xd',:,0,'q',2]) for k in sol.keys()]  # list of max of q[2]
    minheight = [min(sol[k][0]['Xd',:,0,'q',2]) for k in sol.keys()]   # list of min of q[2]
    plt.figure()
    plt.plot(maxheight)
    plt.plot(minheight)
    plt.ylabel('operation heights')
    plt.xlabel('time step (3-hour intervals)')
    plt.legend(['maxheight', 'minheight'])
    plt.title('AWE power output')
    plt.grid('on')

def getWTpower(wind, diam, maxcap, cp = 0.5):
    "input windarray, diameter, maxcap in MW"
    power = []
    for w in wind:
        if w <= 3 or w >= 22:
            power.append(0.)
        else:
            P = 0.5 * 1.225 * cp * np.pi/4 * diam**2 * w**3
            if P*1e-6>=2:
                power.append(maxcap)
            else:
                power.append(P*1e-6)
    return power

    return maxheight, minheight

# --- extract and save the solution
Pcap = '0.666MW'
sol, _, _ = extract_solution(directory, prefix, Pcap)
Pcap = '0.6666MW'
sol, _, _ = extract_solution(directory, prefix, Pcap)
Pcap = '2MW'
sol, _, _ = extract_solution(directory, prefix, Pcap)

maxheight, minheight = operation_altitudes(sol);
sys.exit()


# --- Load the solution
normwind = np.load(directory + prefix+'_features.npy')
heights =  np.load(directory + prefix+'_heights.npy')
timepoints   = normwind.shape[1]
size1 = '0.666MW'; size2 = '0.6666MW'; size3 = '2MW'
kitesizes =  [size1, size2, size3]
plt.ion()
# --- get WT
WTpower       = getWTpower(normwind[-4,:,0],95,2,0.48)

# --- plot load duration and chronical power output
AWEpower = dict()
for cap in kitesizes:
    try:
        with open(directory+'power_'+cap, 'r') as f:
            (AWEpower[cap]) = pickle.load(f)
        plt.figure('load duration curves')
        plt.plot(load_duration_curve(AWEpower[cap]))
        plt.figure('power output')
        plt.plot(AWEpower[cap])
    except:
        print 'Size ', cap , ' not yet available'
legname = [size1, size3]
legname.append('2MW WT')
plt.figure('load duration curves')
plt.plot(load_duration_curve(WTpower))
plt.ylabel('Power in MW')
plt.xlabel('time step (3-hour intervals)')
plt.title('load duration curve')
plt.legend(legname)

plt.figure('power output')
plt.plot(WTpower)
plt.ylabel('Power in MW')
plt.xlabel('time step (3-hour intervals)')
plt.legend(legname)
plt.title('generation')

# --- multiple kites of the < 2.2MW capacity wings
plt.figure('load duration curves ADDED')
plt.plot(load_duration_curve(np.multiply(3,AWEpower[size1])))
plt.plot(load_duration_curve(np.multiply(3,AWEpower[size2])))
plt.plot(load_duration_curve(AWEpower[size3]))
plt.plot(load_duration_curve(WTpower))
plt.ylabel('Power in MW')
plt.xlabel('time step (3-hour intervals)')
plt.title('load duration curve with multiple wings')
plt.legend(['3 x ' + size1, '3 x ' + size1, size3, '2MW WT'])
plt.show()

# # --- zoom in a time period
# fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
# # ax1.plot(np.multiply(3,AWEpower[size1][2100:2300]))
# ax1.plot(np.multiply(2,AWEpower[size2][2100:2300]))
# ax1.plot(AWEpower[size3][2100:2300])
# ax1.plot(WTpower[2100:2300],'--')
# ax2.plot(AWEpower[size1][2100:2300])
# ax2.plot(AWEpower[size1][2100:2300])
# ax2.plot(AWEpower[size3][2100:2300])
# ax2.plot(WTpower[2100:2300], '--')
# fig.tight_layout()
# ax1.set_ylabel('Power in MW')
# ax2.set_xlabel('time step (3-hour intervals)')
# ax1.set_title('power generation zoom')
# ax1.legend(['3 x ' +size1, '2 x ' + size2, size3 , '2MW WT'])
# ax2.legend([size1, size2, size3, '2MW WT'])
# plt.show()

## **kvargs
# ---  Calculate Downtimes and annual generation
for cap in kitesizes:
    AWEpower[cap+'_downtime'] =0
    AWEpower[cap+'_totgen'] =0
    for k in AWEpower[cap]:
        if k <= 0: AWEpower[cap+'_downtime'] += 3
        else: AWEpower[cap+'_totgen'] += 3*k
    print cap+'_downtime', AWEpower[cap+'_downtime'], ' h'
    print cap+'_totgen', AWEpower[cap+'_totgen'], ' MWh in total'
WT_downtime = 0; WT_totgen = 0
for k in WTpower:
    if k<= 0: WT_downtime +=3
    else: WT_totgen += 3 * k
print 'WT_downtime', WT_downtime,' h'
print 'WT_totgen', WT_totgen , ' MWh in total'

# plt.figure()
# y = []; y2 = []; y3 = [];
# for cap in kitesizes:
#     y = np.append(y, AWEpower[cap+'_downtime'])
#     y2= np.append(y2, AWEpower[cap+'_totgen'] )
# y3= np.append(y3, np.multiply(AWEpower[size1+'_totgen'], 3) ); y3= np.append(y3, np.multiply(AWEpower[size2+'_totgen'], 2) ); y3= np.append(y3, AWEpower[size3+'_totgen'] )
# y3 = np.append(y3, WT_totgen)
# y = np.append(y, WT_downtime)
# y2 = np.append(y2, WT_totgen)
# x = np.array(range(len(legname)))
# width = 0.2
# plt.subplot(1,2,1)
# plt.bar(x, y, width, color="red"); plt.ylabel('downtime [h]')
# myticks = legname
# plt.xticks(x,myticks)
# plt.subplot(1,2,2)
# plt.bar(x-0.1, y2, width, color="blue");
# plt.bar(x+0.1, y3, width, color="magenta"); plt.ylabel('total generation in MWh ')
# plt.legend(['single systems' , size1 +' x 3, ' + size2 + ' x 2, ' + size3 + ' x 1'], loc = 4)
# myticks = legname
# plt.xticks(x,myticks)
# plt.show()

# --- WIND SPEEDS ----
plt.ion()
plt.figure()
plt.plot(load_duration_curve(normwind[-4,:,0]))
plt.plot(load_duration_curve(normwind[-5,:,0]))
plt.plot(load_duration_curve(normwind[-6,:,0]))
plt.ylabel('wind in [m/s]')
plt.xlabel('time step (3-hour intervals)')
plt.legend(['~100m','~220m', '~340m'])
plt.title('wind speeds sorted')
plt.grid('on')
plt.ion()

# --- if you want to plot everything together!!!
# import plot_Clusters

#SAVE ALL TO A SINGLE PDF
# import matplotlib.backends.backend_pdf
# pdf = matplotlib.backends.backend_pdf.PdfPages(directory + 'kitecompare_Got16.pdf')
# for fig in xrange(1, plt.figure().number): ## will open an empty extra figure :(
#     pdf.savefig( fig )
# pdf.close()


# #
# plt.ion()
# plt.figure()
# plt.plot(load_duration_curve(normwind[-4,:,0]))
# plt.plot(load_duration_curve(normwind[-5,:,0]))
# plt.plot(load_duration_curve(normwind[-6,:,0]))
# plt.ylabel('wind in [m/s]')
# plt.xlabel('time step (3-hour intervals)')
# plt.legend(['~100m','~220m', '~340m'])
# plt.title('wind speeds sorted')
# plt.grid('on')
# plt.ion()

# plt.ylabel('wind in [m/s]')
# plt.xlabel('time step (3-hour intervals)')
# plt.legend(['~100m','~220m', '~340m'])
# plt.title('wind speeds sorted')
# plt.grid('on')
# plt.ion()
# plt.ion()
# plt.figure()
# plt.plot((normwind[-4,:,0]))
# plt.plot((normwind[-5,:,0]))
# plt.plot((normwind[-6,:,0]))
#
plt.ylabel('wind in [m/s]')
plt.xlabel('time step (3-hour intervals)')
plt.legend(['~100m','~220m', '~340m'])
plt.title('wind speeds sorted')
plt.grid('on')
plt.ion()
