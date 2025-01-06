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
    parser.add_argument('--center_out', type = str, help='', default='./output/main_center.csv')
    parser.add_argument('--label_out', type = str,help='', default='./output/main_label.csv')
    parser.add_argument('--device', type = str, help='', default='cpu')
    parser.add_argument('--eps', type = float, help='')
    parser.add_argument('--minpts', type = int, help='')
    parser.add_argument('--core_out', type=str, default='./output/main_core.csv')
    parser.add_argument('--rad_out', type=str, default='./output/main_rad.csv')
    
    args = parser.parse_args()
    
    whole_data = np.loadtxt(args.data, dtype=float, delimiter=',')
    print('data size = ',whole_data.shape)
    dl = DataLoader(whole_data, args.buffer_size, shuffle=False)
    
    print('Data load completed.')
    start = time()
    
    ## First pass
    print('The first pass:')
    # print(args.eps, args.rho / 2)
    rgda = RGDA(args.eps * args.rho / 2)
    for buffer in tqdm.tqdm(dl):
    # for buffer in dl:
        rgda.Scan(buffer)
    # center = ScalingAlg.StreamingKCenter(dl, args.k, args.alpha, args.m, args.device)
    center = rgda.center
    n_center = center.size()[0]
    
    end = time()
    print('k center completed')
    print('|E|={}'.format(center.size()[0]))
    print('The first pass takes {} s'.format(end - start))
    
    second_start = time()
    np.savetxt(args.center_out, center.detach().cpu().numpy(), fmt='%.2f', delimiter=',')
    
    ## Second Pass
    center_counter = torch.zeros([n_center], dtype=int)
    center_radius = torch.zeros([n_center], dtype=float)
    support_point = [[] for i in range(n_center)]
    center_is_core = torch.zeros_like(center_counter)
    
    print('The second pass')
    for buffer in tqdm.tqdm(dl):
        buffer = buffer.to(args.device)
        data_center_dis = torch.cdist(buffer, center)
        nearest_center_ind = data_center_dis.min(1)[1]
        center_counter = center_counter + (data_center_dis <= torch.ones_like(data_center_dis) * args.eps).sum(0)
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
    if max_radius > args.eps * args.rho / 2:
        print('Warning: center_radius is too large, may fail! {} > {}'.format(max_radius, args.eps * args.rho / 2))
        # sys.exit()
    print('The second pass takes {} s'.format(time()-second_start))
    
    third_start = time()
    # memory = torch.tensor(support_point, dtype=float, device=args.device)
    memory = torch.zeros([center.size()[0] * args.minpts, center.size()[1]], dtype=float, device=args.device)
    i = 0
    for l in support_point:
        if len(l) > 0:
            i -= 1
        for d in l:
            memory[i] = d
            i += 1
    memory = memory[0:i]
    print('|M| = ', i)
    # print('the total memory usage / n = ', (i + center.size()[0]))
    print('the total memory usage / n = ', min((i + center.size()[0])/whole_data.shape[0], 1))
    # memory_counter = torch.zeros([memory.size()[0]], dtype=int)

    ## Third Pass
    print('The third pass')
    for buffer in tqdm.tqdm(dl):
        buffer = buffer.to(args.device)
        data_memory_dis = torch.cdist(buffer, memory)
        memory_counter = (data_memory_dis <= args.eps).sum(0).long()
    
    print('The third pass takes {} s'.format(time()-third_start))
    core_merge_start = time()
    
    center_core = center[torch.nonzero(center_counter > args.minpts).view(-1)]
    memory_core = memory[torch.nonzero(memory_counter > args.minpts).view(-1)]
    core = torch.cat([center_core, memory_core], dim = 0)
    print('|S_*| = ', core.size()[0])
    core_radius = torch.cat([center_radius[torch.nonzero(center_counter > args.minpts).view(-1)], torch.zeros([memory_core.size()[0]], dtype = int)], dim = 0).to(args.device)
    core_label = torch.zeros([core.size()[0]], dtype=int, device=args.device)
    core_core_dis = torch.cdist(core, core, p=2)
    current_label = 1
    # print('Merge inside S_*:')
    for i in range(core.size()[0]):
        if core_label[i] >= 0:
            q = deque()
            q.append(i)
            # print('hello')
            while len(q) != 0:
                # print(len(q))
                p = q.pop()
                core_label[p] = current_label
                p_center_dis = core_core_dis[p]
                neighbor = torch.nonzero(p_center_dis <= core_radius[p] + core_radius + args.eps).view(-1)
                for n in neighbor:
                    if core_label[n] == 0:
                        q.append(n)
                        core_label[n] = -2
            current_label += 1
    print('Number of cluster is', current_label)
    # print('|S_*|=', core.size()[0])
    print('core merge takes {} s'.format(time()-core_merge_start))
    np.savetxt(args.core_out, core.detach().cpu().numpy(), fmt = '%.4f', delimiter = ',')
    np.savetxt(args.rad_out, core_radius.detach().cpu().numpy(), fmt = '%.4f', delimiter = '\n')
    
    fourth_pass_start = time()
    with open(args.label_out, "w+") as f:
        for data in tqdm.tqdm(dl):
            data = data.to(args.device)
            data_core_dis = torch.cdist(data, core)
            nearest_core_dis, nearest_core_ind = data_core_dis.min(1)
            logic_vec = (nearest_core_dis < args.eps + core_radius[nearest_core_ind])
            data_label = logic_vec * (core_label[nearest_core_ind]) + torch.logical_not(logic_vec) * torch.zeros_like(core_label[nearest_core_ind])
            
            for d in data_label:
                f.write(str(int(d)) + '\n')
    print('clustering completed with {} s'.format(time() - start))