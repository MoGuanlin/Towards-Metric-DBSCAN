# Towards Metric DBSCAN: Exact, Approximate, and Streaming Algorithms

This repository is the official implementation of **Towards Metric DBSCAN: Exact, Approximate, and
Streaming Algorithms**.

The Faster Exact Metric DBSCAN Algorithms, and the Metric $\rho$-approximate DBSCAN Via Core Points Summary Algorithms are implemented by c++.

The Streaming $\rho$-approximate DBSCAN Algorithm( Our streaming algorithm ) is implemented by python.

## Exact and approximate algorithms

### Requirements


You can directly run our algorithms implemented by c++ on Linux platform.


### Main Organization of the Code


* **datasets:** This folder contains the clustering datasets we used in the paper. The subdirectory 'GT' contains the datasets only for the baseline 'GT_Exact' and 'GT_Approx' because of the special IO implementation in their source code. 
Please note that if we were to upload all the large datasets we used, it would exceed 100GB of space, which is not feasible for our cloud storage platform. Therefore, we have not uploaded some of the extremely large datasets used in the experiments. Nevertheless, we have provided links to all the datasets in the main text, and you can download them from there.


* **baselines:** This folder contains some open-source code for baselines along with their README files.

* **results:** this folder is used to store the output of the algorithms.
* **src:** this folder contains the C++ code of our algorithms. Besides, folder src/StreamingAlg contains the python code of our streaming algorithms.
* **header:** this folder contains the C++ headers of our algorithms.

### How to Run



**run** file Linux/Algorithms in the following way:

```{Shell}
./linux/Algorithms --alg {algorithm number} --dim {dimension of dataset} --n {number of data points} --minpts {parameter MinPts} --epsilon {parameter epsilon} --rho {approximate ratio rho} --input {dataset file path} --output_data {result output path} --output_center {deprecated}
```

The 'algorithm number' argument can be:

* **1:** vanilla exact DBSCAN 
* **2:** our exact DBSCAN via Radius-Guided Gonzalez
* **3:** our exact DBSCAN via Cover Tree
* **4:** our exact DBSCAN via Radius-Guided Gonzalez for edit distance
* **5:** our approximate DBSCAN via Radius-Guided Gonzalez
* **6:** our approximate DBSCAN via Randomized Gonzalez
* **7:** our approximate DBSCAN via Radius-Guided Gonzalez for edit distance
* **8:** our approximate DBSCAN via Cover Tree
* **9:** vanilla exact DBSCAN for edit distance
* **10:** vanilla exact DBSCAN for angular distance
* **11:** our exact DBSCAN for angular distance
* **12:** our approximate DBSCAN for angular distance



Their are some optional arguments should be used if you run the algorithm 3, 6 or 8:

For Randomized gonzalez-based algorithm:

* **--init_num {initial center number}** 
* **--farthest_num {number of farthest points chosen in every loop}** 
* **--select_num {number of farthest points added to E in every loop}** 

For covertree-based algorithms:
* **--max_layer {max possible covertree layer}** 
* **--terminate_layer {terminate layer when searching NN}** 
  

For example, you can run our approximate method via Radius-Guided Gonzalez by the following command:

```{Shell}
./src/Algorithms --alg 5 --dim 2 --n 4811 --minpts 4 --epsilon 0.03 --rho 0.5 --input ./datasets/banana --output_data ./output_result.txt
```

**compile**

If necessary, you can recompile the executable files in the following way:

```{Shell}
g++ -O3 -std=c++11 -o ./src/Algorithms ./src/*.cpp
```

### Results of Algorithm 1~12

When the user runs Algorithm 1~12, the program will output the following format to the output file specified by the user: 

Each number indicating the clustering result (cluster ID) of one data point.

In addition, the program will output some statistical results to the screen.

For example, if the user enters the following command:

```{Shell}
./src/Algorithms --alg 5 --dim 2 --n 4811 --minpts 4 --epsilon 0.03 --rho 0.5 --input ./datasets/banana --output_data ./output_result.txt
```


and the following content will be output to the screen:

```shell
17 args received.
Running algorithm: 5
load data complete!
*******************K-center procedure*******************
101 loops solved.
max_r = 0.0417732; r = 0.015
201 loops solved.
max_r = 0.0272029; r = 0.015
301 loops solved.
max_r = 0.0206155; r = 0.015
401 loops solved.
max_r = 0.0170294; r = 0.015
Kcenter complete! Centers number: 491
Kcenter time: 24.0225 ms
********************************************************
*******************Calculating neighbors*******************
Calculate Ap complete!
Calculate neighbors time: 0.957504 ms
********************************************************
*******************Mark core procedure*******************
dense count: 360
sparse count: 131
Mark core complete!
Mark core time: 0.251886 ms
********************************************************
*******************Cluster core procedure*******************
Max cluster number is : 2
Cluster core complete!
Cluster core time: 0.347674 ms
********************************************************
*******************Cluster border procedure*******************
Cluster border complete!
Cluster border time: 0.018009 ms
********************************************************
Complete!
Total time: 25.6377 ms
Core point count: 4803
Border point count: 7
Outlier count: 1
```



## Streaming algorithm

Our streaming algorithm is implemented by python.

Before running the program, please ensure that the following libraries are installed in your Python environment.

* PyTorch 2 (cuda is not necessary)
* numpy
* tqdm

To run the example program, please follow these steps:
1. Change the current directory to `src/StreamingAlg`.
2. Enter the command `make moons` to run the program with the `Moons` dataset.
3. The `ARI` and `AMI` indices will be displayed on the screen.
4. If your environment does not support the `make` command, you can open the `Makefile` file with a text editor and manually enter the commands from there into the command line, for example, if you want to run Our streaming algorithm on `Moons` dataset, you need to enter:
    ```{shell}
    moons:
	python ./DBSCAN.py \
		--data=../../datasets/moons_data.csv \
		--buffer_size=100 \
		--rho=0.5 \
		--eps=2 \
		--minpts=10 \
		--label_out=./output/moons_label.csv
	python ./Analysis.py \
		--y=../../datasets/moons_label.csv \
		--label=./output/moons_label.csv
    ```
    

This folder contains five files: 
1. `main.py`: our main code for streaming approximate DBSCAN.
1. `DBSCAN.py`: Another version for streaming approximate DBSCAN.
2. `Analysis.py`: calculate ARI and AMI.
3. `PlotCenter.py`: plot the center of balls.
4. `PlotCore.py`: plot the core points.
5. `PlotLabel.py`: plot the label of points.

### How to run

First, we need to change the current directory to `streamingAlg`.

Then we can run Our streaming algorithm by entering the following command:
```
make moons
```
Then the following result will be displayed in screen.

```
data size =  (9948, 2)
Data load completed.
The first pass: 100%|██████████████████████████████| 78/78 [00:01<00:00, 68.75it/s]
The second pass: 100%|████████████████████████████| 78/78 [00:00<00:00, 110.55it/s]
The third pass: 100%|█████████████████████████████| 78/78 [00:00<00:00, 833.19it/s]
Merge inside S_*: 100%|████████████████████████| 220/220 [00:00<00:00, 7112.44it/s]
|E|= 407
|M| =  1193
|S_*| =  220
Number of clusters  3
the total memory usage / n =  0.16083634901487737
clustering completed with 2.1617186069488525 s
ARI =  0.9675238354642742
AMI =  0.9337253148748574
```
Besides, the program will store the label of all points on src/Algorithm4/output/moons_label.csv

If you want to change data set, please enter the corresponding command in this table.

|    Dataset    | ARI  | AMI  |      Command       |
| :-----------: | :--: | :--: | :----------------: |
|     Moons     | 0.94 | 0.90 |     make moons     |
|    Cancer     | 0.94 | 0.92 |    make cancer     |
|  Arrhythmia   | 0.57 | 0.73 |  make arrhythmia   |
|    Biodeg     | 0.22 | 0.15 |    make biodeg     |
|     MNIST     | 0.70 | 0.52 |     make mnist     |
| Fashion MNIST | 0.47 | 0.71 | make fashion_mnist |
|    USPS HW    | 0.57 | 0.53 |    make uspshw     |
|    Spotify    | 0.02 | 0.14 |    make spotify    |

### For the Million-scale dataset Spotify

Owing to the limitation of the size of submitted files, we can only provide partial dataset of Spotify. The whole dataset needs the size of over 100GB. So we select only 1% of the whole dataset. The data is stored in `/datasets/spotify_session/`.

You can run our streaming algorithm on spotify dataset using the following command:

```shell
make spotify
```

The expected output is 

```shell
python ./DBSCAN4session_faster_4.py \
        --data=./datasets/spotify_session/ \
        --buffer_size=1 \
        --size=2 \
        --rho=0.5 \
        --eps=30.0 \
        --minpts=1 \
        --label_out=./output/session_overall_2_label.csv
Data load completed.
The first pass:
  0%|                                                                                                               | 0/2 [00:00<?, ?it/s]center.shape =  torch.Size([463, 14])
 50%|███████████████████████████████████████████████████                                                   | 1/2 [01:48<01:48, 108.12s/it]center.shape =  torch.Size([610, 14])
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 2/2 [02:29<00:00, 74.67s/it]
k center completed
|E|=610
The first pass takes 149.34458112716675 s
The second pass takes 0.015554666519165039 s
The third pass
The third pass takes 4.0531158447265625e-06 s
|S_*| =  610
Number of cluster is 611
core merge takes 6.7200376987457275 s
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 2/2 [00:38<00:00, 19.47s/it]
clustering completed with 195.0722417831421 s
python ./gen_label.py --data=./datasets/spotify_session/ --size=2 --label_out=./datasets/spotify_session/session_overall_2_label.csv
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 2/2 [00:00<00:00, 71.77it/s]
python ./Analysis.py \
        --y=./datasets/spotify_session/session_overall_2_label.csv \
        --label=./output/session_overall_2_label.csv
ARI =  0.022848635578129835
AMI =  0.14336915115308033
```

### How to run baselines

We also provide the code of baseline algorithms. We implement them using R language. All codes are stored in `/StreamCluster/`.

We need four library, which are `stream`, `magrittr`,`proxy`,`aricode`. `stream` library is the most important. You can refer to `https://cran.r-project.org/web/packages/stream/index.html` to check how to install.

By using R language with `stream` library, we implement four baseline streaming methods mentioned in our paper, which are `DBStream`, `D-Stream`, `evoStream`, `BICO`.

You can run the R code file with the dataset name. For example, if you want to run baselines on `Moons` dataset, you need to run `Moons.r`. The code will call four baselines one by one and output their ARI and AMI.

We have chose the best parameters of every baseline in every dataset. You can also try other parameter by modifying the source code.