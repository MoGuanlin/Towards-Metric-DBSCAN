import numpy as np
import argparse
from sklearn import metrics

parser = argparse.ArgumentParser()
parser.add_argument('--label', type=str,help='data file path')
parser.add_argument('--y', type=str,help='data file path')

args = parser.parse_args()

labels = np.loadtxt(args.label)
y = np.loadtxt(args.y)

print('ARI = ', metrics.adjusted_rand_score(y, labels))
print('AMI = ', metrics.adjusted_mutual_info_score(y, labels))
