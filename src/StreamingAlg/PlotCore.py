import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mp
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--data', type=str,help='data file path')
parser.add_argument('--core', type=str,help='center file path')
parser.add_argument('--rad', type=str,help='center file path')
parser.add_argument('--out', type=str,help='output file path')
parser.add_argument('--data_delimiter', type=str, default=',')
parser.add_argument('--core_delimiter', type=str, default=',')

args = parser.parse_args()

data = np.loadtxt(args.data, delimiter=args.data_delimiter)
centers = np.loadtxt(args.core, delimiter=args.core_delimiter)
rad = np.loadtxt(args.rad, dtype=float)

fig, ax = plt.subplots()

ax.scatter(data[:, 0], data[:, 1], c = 'g')
ax.scatter(centers[:, 0], centers[:, 1], c = 'r')
for i in range(len(rad)):
    a = mp.Circle(centers[i], radius = max(rad[i], 0.001), color = 'b', alpha = 0.5)
    ax.add_artist(a)

# plt.show()
plt.savefig(args.out)