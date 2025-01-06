import sys
from DBSCANPP import DBSCANPP
import numpy as np
import time



""" # Mixture of three multivariate Gaussians
cov = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

X = np.concatenate([np.random.multivariate_normal([5, 5, 5], cov, 1000),
          np.random.multivariate_normal([0, 0, 0], cov, 1000),
          np.random.multivariate_normal([-5, -5, -5], cov, 1000)])
y = np.concatenate([np.full(1000, 0), np.full(1000, 1), np.full(1000, 2)])

# Declare a DBSCAN++ model with tuning hyperparameters
dbscanpp = DBSCANPP(p=0.1, eps_density=5.0, eps_clustering=5.0, minPts=10)
y_hat = dbscanpp.fit_predict(X, init="k-centers")

# Score the clustering
from sklearn.metrics.cluster import adjusted_rand_score, adjusted_mutual_info_score
print("Adj. Rand Index Score: %f." % adjusted_rand_score(y_hat, y))
print("Adj. Mutual Info Score: %f." % adjusted_mutual_info_score(y_hat, y)) """

data = np.loadtxt('{0}'.format(sys.argv[1]), delimiter=' ')
print(data.shape)

np.set_printoptions(threshold=np.inf)
start_time = time.time()
uspspp = DBSCANPP(p=float(sys.argv[2]), eps_density=float(sys.argv[3]), eps_clustering=float(sys.argv[4]), minPts=int(sys.argv[5]))
usps_y = uspspp.fit_predict(data, init="k-centers")
end_time = time.time()
#print(usps_y)
print("time cost: " , end_time-start_time , " s. ")

print(usps_y.shape)
# python test.py