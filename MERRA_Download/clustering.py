#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Clusters the .npy files of wind main and deviation to N clusters.
1. >>merra_download.py<<  : Download wind Data
2. >>merra2npy.py<<       : read files and save as .npy in order to have the needed format for optimisation
3. >>clustering.py<<      : cluster .npy files to get different clusters
this code is executed by merra_download

Created on Mon Nov 20 14:45:55 2017

@author: vive
"""

import numpy as np
import pylab
from sklearn.cluster import KMeans
import matplotlib.cm as cm

# directory = '/Users/elenama/Documents/Git/Pamplona16/'
# prefix ='-1.65_42.81_20160101_20161231_'

def cluster(directory,prefix):
    path = directory+prefix
    timeseries_height = np.load(path+'_heights.npy')
    timeseries_wind   = np.load(path+'_winds.npy')
    timeseries_features = np.load(path+'_features.npy')
    ndays = int(timeseries_wind.shape[1]/8.)
    pylab.figure()
    # --- Plotting raw wind data and extracted features
    pylab.plot(timeseries_height.T[0:8*ndays],'r-')
    #pylab.yscale('log')
    pylab.xlabel('time')
    pylab.ylabel('height (m)')
    pylab.title('Heights of the measurements')
    pylab.show()

    # Feature 1: Absolute wind speed
    pylab.plot(timeseries_features[:,0:8*ndays,0].T,'b-')
    pylab.show()

    # Feature 2: Deviation from angular mean
    pylab.plot(timeseries_features[:,0:8*ndays,1].T,'b-')
    pylab.show()



    # --- Clustering
    nAltitudeLevels = timeseries_features.shape[0]

    # Preprocess data
    clustering_data = np.zeros((ndays*8,2*nAltitudeLevels),dtype='float32')
    for i in range(ndays*8):
        absolutes = timeseries_features[:,i,0]
        delta_radians = timeseries_features[:,i,1]
        clustering_data[i,:] = np.concatenate((absolutes,delta_radians))


    # K-means clustering
    n_clusters = 7
    kmeans = KMeans(n_clusters=n_clusters).fit(clustering_data)
    #kmeans.labels_
    #kmeans.cluster_centers_
    np.save(path+'_kmeans_cluster_labels',kmeans.labels_)
    np.save(path+'_kmeans_cluster_centers',kmeans.cluster_centers_)
    np.save(path+'_clustering_data',clustering_data)



    # --- Plot clustering based on different clusters
    kmeans_labels          = np.load(path+'_kmeans_cluster_labels.npy')
    kmeans_cluster_centers = np.load(path+'_kmeans_cluster_centers.npy')
    clustering_data        = np.load(path+'_clustering_data.npy')
    n_clusters              = kmeans_cluster_centers.shape[0]
    # Assuming wind on altitude levels make up half of the features
    nAltitudeLevels = int(clustering_data.shape[1]/2)
    # https://pythonspot.com/en/matplotlib-bar-chart/

    performance = np.histogram(kmeans_labels,bins=n_clusters)[0]
    y_pos = np.arange(n_clusters)
    pylab.bar(y_pos, performance, align='center', alpha=0.5)
    pylab.xticks(y_pos, y_pos)
    pylab.ylabel('Wind profiles')
    pylab.xlabel('Cluster name')
    pylab.title('Distribution of clusters')
    pylab.show()

    pylab.figure()
    # Plot cluster centers, with respect to absolute wind profiles
    # Plotting: reverse/flip in order of first point=lowest altitude
    cluster_centers_wind_sorted = np.flip(kmeans_cluster_centers[:,0:nAltitudeLevels],1)
    X = np.zeros((n_clusters,nAltitudeLevels))
    Y = np.zeros((n_clusters,nAltitudeLevels))
    for i in range(n_clusters):
        Y[i] = np.arange(nAltitudeLevels)
        X[i] = cluster_centers_wind_sorted[i]
    pylab.plot(X.T,Y.T, 'k')





    # Naive way of plotting: Going through all profiles several times
    # Should only need to pass them once, can be easily improved
    colors = cm.rainbow(np.linspace(0, 1, n_clusters))

    for i in range(n_clusters):
        for j in range(len(kmeans_labels)):
            if kmeans_labels[j] == i:
                #print("wind profile", j, "is in cluster", i)
                # Plotting: reverse/flip in order of first point=lowest altitude
                features = clustering_data[j]
                windprofile_sorted = np.flip(features[0:nAltitudeLevels],axis=0)
                pylab.plot(windprofile_sorted,np.arange(nAltitudeLevels),color=colors[i],alpha=0.2)
    pylab.xlabel('absolute wind speed')
    pylab.ylabel('pressure level')
    pylab.xlim([0,42])
    pylab.title(prefix+'\n'+str(len(kmeans_labels))+' profiles, '+str(n_clusters)+' clusters')
    pylab.savefig(path+'_clusters.pdf', bbox_inches='tight',pad_inches=0.2)

def getClusterwindpairs(directory,prefix):
    path = directory+prefix
    featuresdata = np.load(directory + prefix+'_features.npy')
    cluster_labels  = np.load(directory + prefix+'_kmeans_cluster_labels.npy')
    meanclu = {}
    for c in range(7):
        m = []
        pos = np.where(cluster_labels == c)[0]
        wind = featuresdata[5,pos,0] # sorted at 300 m height (therefore 5)
        for k in wind:
            m.append(np.mean(k))
        meanclu[c] = np.mean(m)
    file = open(directory+'ClusterWindPairs.txt', 'w')
    file.write('Cluster: mean wind speed \n\n')
    file.write(str(meanclu))
    file.close()

# --- EXECUTE
# cluster(directory,prefix)




#plt.plot(cluster_centers_wind.T)
#for cluster_i in range(len(cluster_centers_wind)):
#    cluster = cluster_centers_wind[cluster_i]
#    #indexed_cluster = np.zeros(nAltitudeLevels,dtype='2float32')
#    indexed_cluster = np.zeros(nAltitudeLevels)
#    level = 0
#    for windspeed in cluster:
#        #indexed_cluster[level] = (level,windspeed)
#        indexed_cluster[level] = windspeed
#        level = level + 1


#cluster_labels = np.load('jan2017_kmeans_cluster_labels.npy')
#for i in range(8):
#    for j in range(len(cluster_labels)):
#        if cluster_labels[j] == i:
#            print("wind profile", j, "is in cluster",i)
