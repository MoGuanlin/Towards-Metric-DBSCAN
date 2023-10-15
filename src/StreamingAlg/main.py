import torch
# import ScalingAlg
import argparse
from torch.utils.data import DataLoader
import numpy as np
from time import time
from collections import deque
import tqdm

class RGDA:
    ## Radius guided doubling algorithm
    def __init__(self, r, device = 'cpu') -> None:
        self.center = torch.tensor([], dtype = float).to(device)
        self.r = r
    
    def Scan(self, buffer):
        for data in buffer:
            
            if self.center.size()[0] == 0:
                self.center = torch.cat((self.center, data), dim = 0)[None, :]
                continue
            
            data_center_dis = torch.cdist(data[None, :], self.center, p = 2).view(-1)
            min_dis = data_center_dis.min()
            if(min_dis > self.r):
                self.center = torch.cat((self.center, data[None, :]), dim = 0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type = str,help='data file path')
    parser.add_argument('--buffer_size', type = int)
    parser.add_argument('--rho', type = float, default=2)
    parser.add_argument('--center_out', type = str, help='save the center points', default='./output/main_center.csv')
    parser.add_argument('--label_out', type = str,help='save the label of points', default='./output/main_label.csv')
    parser.add_argument('--device', type = str, help='cpu or cuda', default='cpu')
    parser.add_argument('--r', type = float, help='epsilon in paper', default=100)
    parser.add_argument('--minpts', type = int, help='MinPts in paper')
    parser.add_argument('--core_out', type=str, help='save core points', default='./output/main_core.csv')
    parser.add_argument('--rad_out', type=str, help='save the radius of every ball', default='./output/main_rad.csv')
    
    args = parser.parse_args()
    
    whole_data = np.loadtxt(args.data, dtype=float, delimiter=',')
    X_n = whole_data.shape[0]
    print('data size = ',whole_data.shape)
    dl = DataLoader(whole_data, args.buffer_size, shuffle=False)
    
    print('Data load completed.')
    start = time()
    
    ## First pass
    rgda = RGDA(args.r * args.rho / 2)
    bar = tqdm.tqdm(dl)
    bar.set_description('The first pass')
    for buffer in bar:
        rgda.Scan(buffer)
        
    center = rgda.center
    n_center = center.size()[0]
    
    np.savetxt(args.center_out, center.detach().cpu().numpy(), fmt='%.2f', delimiter=',')
    
    ## Second Pass
    center_counter = torch.zeros([n_center], dtype=int)
    center_radius = torch.zeros([n_center], dtype=float)
    support_point = [[] for i in range(n_center)]
    center_is_core = torch.zeros_like(center_counter)
    bar = tqdm.tqdm(dl)
    bar.set_description('The second pass')
    for buffer in bar:
        # bar.display('The second pass')
        buffer = buffer.to(args.device)
        data_center_dis = torch.cdist(buffer, center)
        nearest_center_ind = data_center_dis.min(1)[1]
        center_counter = center_counter + (data_center_dis <= torch.ones_like(data_center_dis) * args.r).sum(0)
        for i in range(buffer.size()[0]):
            c = nearest_center_ind[i]
            center_radius[c] = max(center_radius[c], data_center_dis[i][c])
            if center_counter[c] < args.minpts and center_is_core[c] == 0:
                support_point[c].append(buffer[i])
            else:
                support_point[c].clear()
            if center_counter[c] >= args.minpts or len(support_point[c]) >= args.minpts:
                center_is_core[c] = 1
    
    max_radius = center_radius.max()

    memory = torch.zeros([center.size()[0] * args.minpts, center.size()[1]], dtype=float, device=args.device)
    i = 0
    for l in support_point:
        for d in l:
            memory[i] = d
            i += 1
    M_size = i
    memory = memory[0:i]

    ## Third Pass
    bar = tqdm.tqdm(dl)
    bar.set_description('The third pass')
    for buffer in bar:
        buffer = buffer.to(args.device)
        data_memory_dis = torch.cdist(buffer, memory)
        memory_counter = (data_memory_dis <= args.r).sum(0).long()
    
    center_core = center[torch.nonzero(center_counter > args.minpts).view(-1)]
    memory_core = memory[torch.nonzero(memory_counter > args.minpts).view(-1)]
    core = torch.cat([center_core, memory_core], dim = 0)
    core_radius = torch.cat([center_radius[torch.nonzero(center_counter > args.minpts).view(-1)], torch.zeros([memory_core.size()[0]], dtype = int)], dim = 0).to(args.device)
    core_label = torch.zeros([core.size()[0]], dtype=int, device=args.device)
    core_core_dis = torch.cdist(core, core, p=2)
    current_label = 1
    bar = tqdm.tqdm(range(core.size()[0]))
    bar.set_description('Merge inside S_*')
    for i in bar:
        if core_label[i] == 0:
            q = deque()
            q.append(i)
            while len(q) != 0:
                p = q.pop()
                core_label[p] = current_label
                p_center_dis = core_core_dis[p]
                neighbor = torch.nonzero(p_center_dis <= core_radius[p] + core_radius + args.r).view(-1)
                for n in neighbor:
                    if core_label[n] == 0:
                        q.append(n)
                        core_label[n] = -2
                        
            current_label += 1
    np.savetxt(args.core_out, core.detach().cpu().numpy(), fmt = '%.4f', delimiter = ',')
    np.savetxt(args.rad_out, core_radius.detach().cpu().numpy(), fmt = '%.4f', delimiter = '\n')
    
    with open(args.label_out, "w+") as f:
        bar = tqdm.tqdm(dl)
        for data in bar:
            data = data.to(args.device)
            data_core_dis = torch.cdist(data, core)
            nearest_core_dis, nearest_core_ind = data_core_dis.min(1)
            logic_vec = (nearest_core_dis < args.r + core_radius[nearest_core_ind])
            data_label = logic_vec * (core_label[nearest_core_ind]) + torch.logical_not(logic_vec) * torch.zeros_like(core_label[nearest_core_ind])
            
            for d in data_label:
                f.write(str(int(d)) + '\n')
    
    print('|E|=', n_center)
    print('|M| = ',M_size)
    print('|S_*| = ', core.size()[0])
    print('Number of clusters ', current_label)
    print('the total memory usage / n = ', min((M_size + n_center)/X_n, 1))
    print('clustering completed with {} s'.format(time() - start))