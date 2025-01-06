
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--data', type=str,help='data file path')
parser.add_argument('--label', type=str,help='label file path')
parser.add_argument('--out', type=str,help='output file path')
args = parser.parse_args()

data = np.loadtxt(args.data, delimiter=',')
labels = np.loadtxt(args.label, dtype=int)

num_labels = labels.max() + 1
data_color = []
for color in range(-1, num_labels):
    new_color = np.nonzero(labels == color)
    data_color.append(new_color)

color_table = {-1:'r', 0:'k', 1:'g', 2:'b', 3:'y', 4:'brown', 5:'cyan', 6:'tan', 7:'tomato', 8:'violet',
               9:'#FFB6C1', 10:'#FFC0CB', 11:'#DC143C', 12:'#FF0F55', 13:'#DB7093',
               14:'#FF69B4', 15:'#FF1493', 16:'#C71585', 17:'#DA70D6', 18:'#D8BFD8',
               19:'#FABCD4', 20:'#FF1502', 21:'#CC1502', 22:'#7FF2E1', 23:'#A8F63D'}

len(data_color)

for color in range(num_labels + 1):
    plt.scatter(data[data_color[color]][:, 0], data[data_color[color]][:, 1], color = color_table[color - 1]) 

plt.savefig(args.out)
