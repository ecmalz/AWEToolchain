# -- Import the OCP
from init_Drag_wind import *
import time
# opt, val_opt = create_initarg(6.5,1,vars_init)
#
# sys.exit()

# Creating the init-files if non-existent.
wind_init = [4, 6, 10 ,12]
free_l = 0 # put 1 if you want to free the tether for arg_init
init = vars_init

Pcap = params['Pnom']
import os.path
for w in wind_init:

    if os.path.isfile('arg_Pcap' + str(w) + '.dat'):
        print '\n\n  arg_Pcap' + str(w) + '.dat is existent. '

    else:
        print '\n\n arg_Pcap'+ str(w) + '.dat has to be created by running the OCP with wind speed ' + str(w) + '. '
        time.sleep(1)
        create_initarg(w,free_l,init)

time.sleep(1)

choice = raw_input('Do you wanna exit now (e) or proceed (p)? ')
if choice == 'e':
    sys.exit()
else:
    pass





print '#################################################'
print 'Start running through different wind clusters ... '
print '#################################################'

cities = {'Pamplona': '42.81_-1.65', 'Esbjerg': '55.50_8.50', 'Goteborg': '57.5_11.875', 'Paris': '48.86_2.35',\
    'laPalma':'28.75_-17.99', 'Galway':'53.27_9.05', 'redsea':'21.73_36.63', 'zermatt':'46.00_7.75'}

# location = 'gothenburg'
case = 'wind'
clu = 7
city = 'laPalma'
directory = '/Users/elenama/Documents/Git/'+city+'16/'
# prefix = '57.5_11.875_20160101_20161231'
prefix = cities[city]+'_20160101_20161231'

# --------------------------------------
# RUN THROUGH DIFFERENT WIND PROFILES
# --------------------------------------
if case=='wind':

    # --- save in correct directory
    import datetime
    # today =  str(datetime.datetime.now().year-2000) + str(datetime.datetime.now().month).zfill(2) +  str(datetime.datetime.now().day).zfill(2)
    today = str(Pcap) + 'MW'
    import os
    file_path = directory + today
    try:
        os.stat(file_path)
    except:
        os.mkdir(file_path)

    import time
    import sys
    start_time       = time.time()
    barrierparameter = 1e-4
    # -- update solver with fixed barrier parameter
    opts['ipopt.mu_init']     = barrierparameter
    opts['ipopt.mu_max']      = barrierparameter
    opts['ipopt.mu_min']      = barrierparameter
    opts['ipopt.mu_target']   = barrierparameter
    opts["ipopt.max_iter"]    = 500  # More will skip that profile
    opts["expand"]            = True
    opts['ipopt.linear_solver'] = 'ma57'
    solver = ca.nlpsol("solver", "ipopt", nlp, opts)
    sys.path.append('../../..Data')

    # winddata        = np.load('../../../Data/'+location+'_20160101_20161231_wind_lowest10.npy')
    # heightsdata     = np.load('../../../Data/'+location+'_20160101_20161231_heights_lowest10.npy')
    # featuresdata    = np.load('../../../Data/'+location+'_20160101_20161231_features_lowest10.npy')
    # cluster_labels  = np.load('../../../Data/'+location+'_20160101_20161231_kmeans_cluster_labels.npy')
    winddata     = np.load(directory + prefix+'_winds.npy')
    heightsdata  = np.load(directory + prefix+'_heights.npy')
    featuresdata = np.load(directory + prefix+'_features.npy')
    cluster_labels  = np.load(directory + prefix+'_kmeans_cluster_labels.npy')
    winddata     = winddata[-10:,:,:]
    heightsdata  = heightsdata[-10:,:]
    featuresdata = featuresdata[-10:,:,:]

    failed          = {}   # Collecting data stamps where the solver failed

    for clu in range(clu):
        pos             = np.where(cluster_labels == clu)[0]  # pos = actual time stamps included in current cluster
        sorted_index    = np.argsort(featuresdata[5,pos,0])   # Sort once more within a cluster - wind at height 300m from low to high, output: sorted index of pos
        sorted_pos      = [pos[k] for k in sorted_index]      # arrange pos with new indexing to sorted_pos

        # --  Main direction of wind = x direction, derivation in y
        xwind = []; ywind = []; heights = []
        for k in sorted_pos:
            x = [w * np.cos(-a) for w,a in featuresdata[:,k,:]]
            y = [w * np.sin(-a) for w,a in featuresdata[:,k,:]]
            # from now cluster-wise sorted defines the new indices
            xwind.append(np.array(x, dtype = float))
            ywind.append(np.array(y, dtype = float))
            heights.append(heightsdata[:,k])

        # --- Choose correct init log profile
        if clu==1 or clu==5:
            init_wind = 4
        elif clu == 3 or clu == 0:
            init_wind = 6
        elif clu ==4 or clu==2:
            init_wind = 10
        if clu ==6:
            init_wind = 12

        # --- open solver data from homotopy (log) solve for windspeed wind0
        with open('arg_Pcap'+ str(init_wind)+ '.dat','r') as f:
            (arg_init, res_init, p_num) = pickle.load(f)

        print '\n\n >>>>>Loop over wind profiles in cluster', clu, '<<<<<<'

        # --- Loop over wind profiles
        h   = range(0,900, 10)
        arg = arg_init
        res = res_init
        arg['lbg']['h','max_power'] =  -Pcap*1e6
        #arg['lbx']['l'] = 400;  arg['ubx']['l'] = 400;
        for k in range(len(sorted_pos)):  # k is NOT rerpesenting time stamp but ind in cluster-wise sorted

            _,tau_x = fcts.smooth_Lagrange_poly(heights[k],xwind[k]) # get actual x poly_parameters
            _,tau_y = fcts.smooth_Lagrange_poly(heights[k],ywind[k]) # get actual y poly_parameters

            p_num['puni','tau_x']    = tau_x        # update parameters in NLP
            p_num['puni','tau_y']    = tau_y        # update parameters in NLP
            p_num['puni','altitude'] = heights[k]

            arg['x0']   = V(res['x'])               # give new initial guess & parameters to solve
            arg['p']    = p_num
            print '####################### SOLVING for Cluster', clu, ' number', k, 'timestep',sorted_pos[k], '#############'
            # --- Solve
            res = solver(**arg)

            # --- Exit if no solution is found
            stats = solver.stats()
            if stats['return_status'] == 'Maximum_Iterations_Exceeded' or stats['return_status'] == 'Restoration_Failed' or stats['return_status'] == 'Solved_To_Acceptable_Level':
                arg = arg_init
                res = res_init
                failed[clu] = sorted_pos[k]
                print 'failed'
                os.system('tput bel && tput bel && tput bel && tput bel && tput bel  ')
                import time; time.sleep(2)
                continue



            # --- Save solution
            opt = V(res['x'])
            # print '###########################################'
            print '\n ############## SOLVED for Cluster', clu, ' number', k, 'timestep',sorted_pos[k], '#############'
            # # m = params['mK'] + 1./3 *   params['tether_density'] * np.mean(opt['Xd',:,:,'ltet'])
            average_power = -opt['Xd',-1,-1,'E']/float(params['ScalePower'])/opt['tf']
            print '\n ---> Average Power = ', average_power, '  Orbit period = ', opt['tf']
            # print 'optimal length', V(res['x'])['l']
            # print '###########################################'

            outputs  = Output(Out_fun(V,P))
            val_init = get_Output(vars_init,p_num)
            val_opt  = get_Output(opt,p_num)
            # print val_opt['speed']
            # print val_opt['windspeed_shear']
            with open(file_path+'/timestep' + str(sorted_pos[k])  + '_'+str(Pcap)+'MW','w') as f:
                pickle.dump((average_power, clu, heights[k], tau_x, tau_y, res),f)
            print '\n ############## timestep ',sorted_pos[k], ' saved !!  #############'


            # --- Plot data approximation
            # import matplotlib.pylab as plt
            # plt.ion()
            # Lagrange_polynomial_x = fcts.Lagrange_poly(heights[k],tau_x)
            # Lagrange_polynomial_y = fcts.Lagrange_poly(heights[k],tau_y)
            # plt.figure('polynomial')
            # plt.subplot(1,2,1)
            # plt.scatter( xwind[k], heights[k])
            # plt.plot([Lagrange_polynomial_x([j])[0] for j in h ], h)
            # plt.ylim([0,800])
            # plt.grid('on')
            # plt.subplot(1,2,2)
            # plt.scatter(ywind[k],heights[k])
            # plt.plot([Lagrange_polynomial_y([j])[0] for j in h ],h)
            # plt.grid('on')
            # plt.ylim([0,800])
            # plt.show()


        print '\n >>>>>> Cluster', clu, 'done! <<<<<< \n'
    end_time = time.time()
    time_taken = end_time - start_time
    print time_taken
    with open(file_path+'/failed','w') as f:
        pickle.dump((failed),f)



if case=='NOcluster_wind':
    location = 'vimmerby'

    # --- save in correct directory
    import datetime
    today =  str(datetime.datetime.now().year-2000) + str(datetime.datetime.now().month) +  str(datetime.datetime.now().day)
    import os
    file_path = "../../../Data/drag_dat/"+location+'12/'+today+'/'
    directory = os.path.dirname(file_path)
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)

    import time
    start_time = time.time()
    barrierparameter =  1e-4
    opts['ipopt.mu_init']   = barrierparameter
    opts['ipopt.mu_max']    = barrierparameter
    opts['ipopt.mu_min']    = barrierparameter
    opts['ipopt.mu_target'] = barrierparameter
    opts["expand"]          = True
    opts['ipopt.linear_solver'] = 'ma57'
    solver = ca.nlpsol("solver", "ipopt", nlp, opts)

    # Get wind data featuresdata[heights, profile/timestamp, [norm,angle deviation]]
    # winddata     = np.load('../../jan2017_wind_lowest10.npy')
    # heightsdata  = np.load('../../jan2017_heights_lowest10.npy')
    # featuresdata = np.load('../../jan2017_features_lowest10.npy')
    # cluster_labels  = np.load('../../jan2017_kmeans_cluster_labels.npy')
    # cluster_labels  = np.load('../../../Data/'+location+'_20160101_20161231_kmeans_cluster_labels.npy')
    # winddata  = np.load('../../../Data/2016_diffcities_lowest15/ts_' + location +'_20160101_20161231_wind_lowest15.npy')
    # featuresdata  = np.load('../../../Data/2016_diffcities_lowest15/ts_' + location +'_20160101_20161231_features_lowest15.npy')
    # heightsdata  = np.load('../../../Data/2016_diffcities_lowest15/ts_' + location +'_20160101_20161231_heights_lowest15.npy')
    # winddata        = np.load('../../../Data/'+location+'_20160101_20161231_wind_lowest10.npy')
    # heightsdata     = np.load('../../../Data/'+location+'_20160101_20161231_heights_lowest10.npy')
    # featuresdata = np.load('../../../Data/'+location+'_20160101_20161231_features_lowest10.npy')

    features = np.load('/Users/elenama/Box Sync/Code/Kitemodel/Data/2016GotVimmerby/ts_'+location+'_20120101_20121231_features_lowest15.npy')
    heightsdata = np.load('/Users/elenama/Box Sync/Code/Kitemodel/Data/2016GotVimmerby/ts_'+location+'_20120101_20121231_heights_lowest15.npy')
    nlevels, ntimesteps,ndata = features.shape
    featuresdata = features[range(5,nlevels),:,:]
    heightsdata = heightsdata[range(5,nlevels),:]

    # --- open solver data from homotopy (log) solve for windspeed wind0
    with open('../../../Data/drag_dat/arg_1MW_' +str(6)+ '.dat','r') as f:
        (arg, res, p_num) = pickle.load(f)
    arg_init = arg
    for k in range(ntimesteps):   # for all the profiles#
    # for k in range(2107, 2108):   # for all the profiles#
    # for k in range(1):   # for all the profiles#

        xwind = [w * np.cos(-a) for w,a in featuresdata[:,k,:]]
        ywind = [w * np.sin(-a) for w,a in featuresdata[:,k,:]]
        xwind = np.array(xwind, dtype = float)
        ywind = np.array(ywind, dtype = float)
        heights = heightsdata[:,k]
        _,tau_x =  fcts.smooth_Lagrange_poly(heights,xwind)
        _,tau_y =  fcts.smooth_Lagrange_poly(heights,ywind)

        p_num['puni','tau_x']    = tau_x        # update parameters in NLP
        p_num['puni','tau_y']    = tau_y        # update parameters in NLP
        p_num['puni','altitude'] = heights

        arg['x0']   = V(res['x'])               # give new initial guess & parameters to solver
        arg['p']    = p_num

        print '#############################'
        # --- Solve
        res = solver(**arg)

        stats = solver.stats()
        if stats['return_status'] == 'Maximum_Iterations_Exceeded' or stats['return_status'] == 'Restoration_Failed':
            arg = arg_init
            res = res_init

        # --- Save solution
        opt = V(res['x'])
        # outputs  = Output(Out_fun([V,P])[0])
        # val_init = get_Output(vars_init,p_num)
        # val_opt  = get_Output(opt,p_num)


        print '###########################################'
        print '\n ############## solved time stamp', k, '#############'
        # m = params['mK'] + 1./3 *   params['tether_density'] * np.mean(opt['Xd',:,:,'ltet'])
        average_power = -opt['Xd',-1,-1,'E']/float(ScalePower)/opt['tf']
        print '\n ---> Average Power = ', -opt['Xd',-1,-1,'E']/float(ScalePower)/opt['tf'], '  Orbit period = ', opt['tf']
        print 'optimal length', V(res['x'])['l']
        print '###########################################'

        if average_power<0:
            arg = arg_init


        with open(directory+'/timestep' + str(k) +'_1MW.dat','w') as f:   # Note: AUTOMATE SAVING
            pickle.dump((average_power, heights, tau_x, tau_y, res),f)


    end_time = time.time()
    time_taken = end_time - start_time
    print time_taken


print 'you are done with the location ' + city + ' and kitesize ' + today + '.  !!! \n\n Lets check your results by running sol_ClusterRun.py. But make sure you load the correct files.'
