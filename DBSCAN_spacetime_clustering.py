import numpy as np
from numpy import matlib as mb
import copy
import numpy.ma as ma
import pandas as pd
import netCDF4 as nc
import scipy
import argparse

from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score as ssc
from sklearn.metrics import davies_bouldin_score as dbsc
from sklearn.metrics import calinski_harabasz_score as chsc



def parse_commandline():
    """Parse the arguments given on the command-line.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year",
                       help="input year of data",
                       default=None)
    args = parser.parse_args()
    return args


def reshape_labs(lab, dshape, vinds, nd):
    labsfull = np.empty(dshape)
    labsfull[:] = np.nan
    labsfull[:, vinds[0], vinds[1]] = np.reshape(lab, (nd, -1))
    return labsfull
        
#%%
    
def make_STdata(data0, lons0, lats0, ss, Lrad, Tday):
    # Stabilizes the distances

    # select spatial data (dim2)
    # data0 = data0[:,:,150:350]
    # lats0 = lats0[:,150:350]
    # lons0 = lons0[:,150:350]

    # # subsample
    data0 = data0[:, ::ss, ::ss]
    lons0 = lons0[::ss, ::ss]
    lats0 = lats0[::ss, ::ss]
    
    lonsav = copy.copy(lons0)
    latsav = copy.copy(lats0)

    mask0 = data0.mask;
    data0 = data0.data;


    lons0 = lons0 - np.min(lons0.flatten())
    lons0 = (np.pi/180)*lons0
    lats0 = (np.pi/180)*lats0

    valinds = np.nonzero(~mask0[0, :, :])
    lons0 = lons0[valinds]
    lats0 = lats0[valinds]
    sz1 = lats0.shape[0]
    
    data0clean = np.reshape(data0[:, valinds[0], valinds[1]].flatten(), (-1, 1))

    dlist = []
    
    numdays = data0.shape[0];

    dlist = (sz1*numdays)*[None]  
    sztemp = 5*round(sz1/100)
    
    for i in np.arange(sz1):
        
        if i%sztemp == 0:
            print("Spatial lagging : "+str(round(5*i/sztemp))+" %")
        
        dlats = lats0[i] - lats0
        val1 = np.multiply(np.sin(lons0[i]),np.sin(lons0)) + np.multiply(np.multiply(np.cos(lons0[i]), np.cos(lons0)), np.cos(dlats))
        val1[val1 > 1] = 1
        mydists = np.arccos(val1)
        mydists[i] = 0 # To numerically stabilize the self distance 
        inds = (mydists <= Lrad)
        
        # print(i)
        
        for j in np.arange(numdays):
            d1 = data0[j, valinds[0], valinds[1]]    
            dlist[j*sz1 + i] = d1[inds]

    dlist2 = []
    tinds = np.kron(np.arange(numdays), np.ones(d1.shape[0]))

    ld = len(dlist)
    ldtemp = 5*round(ld/100)
    for j in np.arange(ld):
        
        if j%ldtemp == 0:
            print("Temporal lagging : "+str(round(5*j/ldtemp))+" %")
        
        inds = np.arange(j - Tday*d1.shape[0], j + Tday*d1.shape[0] + 1, d1.shape[0])
        inds = inds[inds > -1]
        inds = inds[inds < ld]
        inds[inds != j]
        
        dlist2.append(dlist[j])
        
        for k in inds:
            dlist2[j] = np.concatenate((dlist2[j], dlist[k]))
            
        dlist2[j] = np.reshape(dlist2[j].flatten(), (-1, 1))
        
    return data0clean, dlist2, data0.shape, valinds, numdays, lonsav, latsav
  
class K_Means_lagged_mini:
    def __init__(self, k=4, tol=1e-3, maxits=200):
        self.k = k
        self.tol = tol
        self.maxits = maxits
        self.centroids = np.zeros((1, k))

    def fit_predict(self, data, ndpts, bs):

        newcentroids = np.zeros((1, self.k));
        
        batchinds = np.arange(ndpts)
        nbatches  = int(np.floor(ndpts/bs))
        
        if ndpts/bs != np.floor(ndpts/bs):
            nbatches = nbatches + 1
        
        for its in range(self.maxits):
            np.random.shuffle(batchinds)
            
            for batchits in range(nbatches):
                
                mbinds = np.array(batchinds[ np.arange(batchits*bs, np.min(((batchits+1)*bs, ndpts))) ])
                
                batchdata = [data[i] for i in mbinds]

                newcentroids = 0*newcentroids
                classification = self.predict(batchdata, len(batchdata))
                
                for i in range(self.k):
                
                    inds = np.nonzero(classification == i)[0]
                    newcentroids[:, i] = np.mean(np.concatenate([batchdata[j] for j in inds]).flatten())
    
                flag = np.linalg.norm(self.centroids - newcentroids)
                
                print('Iter = {} KM err = {}'.format(its, flag))
                self.centroids = newcentroids
                if flag < self.tol:
                    break
        
        classification = self.predict(data, ndpts)
        return classification

    def predict(self, data, ndpts):
        
        classification = np.zeros(ndpts)
        
        for j in range(ndpts):
            classification[j] = np.argmin(np.mean(np.power(data[j] - self.centroids, 2), 0))
    
        return classification
  
file_path = "/projects/p31970/svdi_drought/data/NLDAS/annual_data/NLDAS_SVDI_"
args = parse_commandline()
year = int(args.year)
f1 = nc.Dataset(file_path + str(year) + ".nc",'r')
data = f1.variables['SVDI'][:]


f2 = nc.Dataset('/projects/p31970/svdi_drought/data/NLDAS_lon_lat.nc','r')
lats = (f2.variables['lat'][:].data)
lons = (f2.variables['lon'][:].data)

lats = np.transpose(mb.repmat(lats, lons.shape[0], 1))
lons = mb.repmat(lons, lats.shape[0], 1)

# Subsampling by 5
ss = 10
# Space clustering around 0.02 radians ~  120 km
Lrad = 0.02
# 7 day time lag on each side
Tday = 7 

### Space-Time Neighborhood creation ###
[dataclean, dataST, datashape, valinds, ndays, lonsp, latsp] = make_STdata(data, lons, lats, ss, Lrad, Tday)


##### K- Means #####
# Number of clusters for space-time kmeans
k = 5

km = MiniBatchKMeans(n_clusters=k)
km.fit(dataclean)
km.cluster_centers_ = np.sort(km.cluster_centers_, 0)
labkm = km.predict(dataclean)
labskmfull = reshape_labs(labkm, datashape, valinds, ndays)

#%%

kmeans_w_lag = K_Means_lagged_mini(k=k, tol=1e-3, maxits=2)
kmeans_w_lag.centroids = np.reshape(km.cluster_centers_, (1, k))

## cluster labels for k-means 
labkmst = kmeans_w_lag.fit_predict(dataST, len(dataST), bs=10000)
# reshape labels to the shape of original data
labskmstfull = reshape_labs(labkmst, datashape, valinds, ndays)
## copy of cluster to manipulate for cluster consolidation 
copied_cluster_ind = labskmstfull.flatten().copy()

# For each cluster, calculate mean and standard deviation
# Combine clusters (by changing their label) if they have means within the 
# same SVDI category
for clusternum in np.unique(labskmstfull):
  cluster_svdi = (reshape_labs(dataclean, datashape, valinds, ndays).flatten())[(labskmstfull.flatten() == clusternum)]
  cluster_mean = np.nanmean(cluster_svdi)
  cluster_sd = np.nanstd(cluster_svdi)
  print('Cluster Number :' + str(clusternum))
  print('Cluster Mean: ' + str(np.nanmean(cluster_svdi)))
  print('Cluster Standard Dev: ' + str(np.nanstd(cluster_svdi)))
  if cluster_mean >= 2:
      # Extreme drought 
    copied_cluster_ind[copied_cluster_ind == clusternum] = k-1
  elif cluster_mean < 2 and cluster_mean >= 1.5:
        # Severe Drought
    copied_cluster_ind[copied_cluster_ind == clusternum] = k-2
  elif cluster_mean < 1.5 and cluster_mean >= 1:
        # Moderate Drought (b/w 1 and 1.5)
    copied_cluster_ind[copied_cluster_ind == clusternum] = k-3
  elif cluster_mean < 1 and cluster_mean >=.5:
        # Mild Drought
    copied_cluster_ind[copied_cluster_ind == clusternum] = k-4
  else:
        # Not dry enough to be considered drought
    copied_cluster_ind[copied_cluster_ind == clusternum] = np.nan 
    
    
# For each drought category, select the cluster associated with that drought category
Ts = np.kron(np.reshape(np.arange(ndays), (-1, 1)), np.ones((np.prod(lonsp.shape), 1)))
Xs = np.kron(np.ones((ndays, 1)), np.reshape(latsp.flatten(), (-1, 1)))
Ys = np.kron(np.ones((ndays, 1)), np.reshape(lonsp.flatten(), (-1, 1)))

pureST0 = np.concatenate((Ts, Xs, Ys), axis=1)
pureST = (pureST0 - np.mean(pureST0, axis=0))/np.std(pureST0, axis=0)

dryST0_extreme = pureST0[(copied_cluster_ind.flatten() == k-1), :]
dryST_extreme = pureST[(copied_cluster_ind.flatten() == k-1), :]

dryST0_severe = pureST0[(copied_cluster_ind.flatten() == k-2), :]
dryST_severe = pureST[(copied_cluster_ind.flatten() == k-2), :]

dryST0_moderate = pureST0[(copied_cluster_ind.flatten() == k-3), :]
dryST_moderate = pureST[(copied_cluster_ind.flatten() == k-3), :]

dryST0_mild = pureST0[(copied_cluster_ind.flatten() == k-4), :]
dryST_mild = pureST[(copied_cluster_ind.flatten() == k-4), :]

#### DBSCAN #####
# select DBSCAN espilon parameter
epsilon_val = .25

extreme_sev_mod_mild = [dryST_extreme,dryST_severe, dryST_moderate,dryST_mild]
extreme0_sev0_mod0_mild0 = [dryST0_extreme,dryST0_severe,dryST0_moderate,dryST0_mild]


clt_extreme_sev_mod_mild = [[],[],[],[]]
labs_extreme_sev_mod_mild = [[],[],[],[]]
i = 0
for drought in extreme_sev_mod_mild:
    num_drought_points = drought.shape[0]
    print(num_drought_points)
    # select only drought clusters with points, if it's empty cluster, there is no
    # drought cluster in that particular category
    if num_drought_points == 0:
        clt_extreme_sev_mod_mild[i] = np.nan
    else:
      # Applies DBSCAN to non-empty clusters
        clt = DBSCAN(eps=epsilon_val).fit(drought)
        clt_extreme_sev_mod_mild[i] = clt.labels_
        labs_extreme_sev_mod_mild[i] = np.ones(clt.labels_.shape, dtype=bool)
    i+=1

# Trims DBSCAN clusters to be above trim value.
for clt_labels_emm,labs_emm,ind in zip(clt_extreme_sev_mod_mild,labs_extreme_sev_mod_mild,np.arange(4)):
    if len(labs_emm) == 0:
        continue
    else:
        labsu = np.unique(clt_labels_emm)
        vols = np.zeros(labsu.shape[0])
        trimmed = 0
        trimvolume = 500
        for i in np.arange(labsu.shape[0]):
            inds2 = (clt_labels_emm == labsu[i])
            if labsu[i] == -1:
                labs_emm[inds2] = False
            else:
                vols[i] = np.sum(inds2)
                if vols[i] < trimvolume:
                    trimmed = trimmed + 1 
                    labs_emm[inds2] = False
                    
cluster_svdi_extreme_sev_mod_mild = [[],[],[],[]]
i=0
# 4 corresponds t extreme drought
# 3 corresponds to severe
# 2 corresponds to moderate
# 1 to mild
for cluster_severity in [4,3,2,1]:
    cluster_svdi = (reshape_labs(dataclean, datashape, valinds, ndays).flatten())[(copied_cluster_ind.flatten() == cluster_severity)]
    print(cluster_severity)
    print('SVDI VALUES SHAPE:' + str(cluster_svdi.shape))
    if cluster_svdi.shape[0] == 0:
        cluster_svdi_extreme_sev_mod_mild[i] = 0
    else:
        cluster_svdi_extreme_sev_mod_mild[i] = cluster_svdi
        print(cluster_svdi.shape)
    i+=1
    
drought_category = ['extreme','severe','moderate', 'mild']
svdi_all = []
day_all = []
lat_all = []
long_all = []
cluster_all = []
drought_severity_all = []
year_all = []
for category,cluster_svdi_emm, dryST0_emm,clt_labels_emm,labs_emm in zip(drought_category,cluster_svdi_extreme_sev_mod_mild, extreme0_sev0_mod0_mild0,clt_extreme_sev_mod_mild,labs_extreme_sev_mod_mild):
    if(np.all(cluster_svdi_emm == 0)):
        continue
    else:
        for cluster_id in np.unique(clt_labels_emm[labs_emm]):
            print('Drought Category: ' + str(category))
            print('Cluster ID: ' + str(cluster_id))
            svdi = cluster_svdi_emm[labs_emm][clt_labels_emm[labs_emm] == cluster_id]
            print('SVDI Length: ' + str(len(svdi)))
            day = (dryST0_emm[labs_emm,0][clt_labels_emm[labs_emm] == cluster_id])
            print('Day Length: ' + str(len(day)))
            lat = (dryST0_emm[labs_emm,1][clt_labels_emm[labs_emm] == cluster_id])
            long = (dryST0_emm[labs_emm,2][clt_labels_emm[labs_emm] == cluster_id])
            cluster_id_val= (np.repeat(np.array(cluster_id),len(day), axis = 0)) 
            drought_severity = (np.repeat(np.array(category),len(day), axis = 0))  
            year_ = np.repeat(np.array(year),len(day), axis = 0) 
            svdi_all.extend(svdi)
            day_all.extend(day)
            lat_all.extend(lat)
            long_all.extend(long)
            cluster_all.extend(cluster_id_val)
            drought_severity_all.extend(drought_severity)
            year_all.extend(year_)
        year_ = np.repeat(np.array(year),len(cluster_id_val), axis = 0) 
        clustering_output = {'day': day_all, 'long': long_all, 'lat': lat_all,'SVDI': svdi_all,'cluster_id': cluster_all, 'drought_severity': drought_severity_all, 'year':year_all} 
        dbscan_output = pd.DataFrame(data=clustering_output) 

file_prefix = 'DBSCAN_output_esmm_' 
dbscan_output.to_csv(file_prefix + str(year)  +".csv")
