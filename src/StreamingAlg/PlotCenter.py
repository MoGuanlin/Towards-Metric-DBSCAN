import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--data', type=str,help='data file path')
parser.add_argument('--center', type=str,help='center file path')
parser.add_argument('--out', type=str,help='output file path')
parser.add_argument('--data_delimiter', type=str, default=',')
parser.add_argument('--center_delimiter', type=str, default=',')

args = parser.parse_args()

data = np.loadtxt(args.data, delimiter=args.data_delimiter)
centers = np.loadtxt(args.center, delimiter=args.center_delimiter)

plt.scatter(data[:, 0], data[:, 1], c = 'g')
plt.scatter(centers[:, 0], centers[:, 1], c = 'r')

# plt.show()
plt.savefig(args.out)