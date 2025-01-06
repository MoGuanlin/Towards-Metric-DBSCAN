
***********************************************************************************
 The DBSCAN program was compiled under Ubuntu 14.04 using g++ with O3 turned on.
*********************************************************************************** 


===================================================================================
 Files
===================================================================================
Unzip the "DBSCAN_Bin.zip" file, which will produce a folder containing:

- a binary file named "DBSCAN",
 
- a sample dataset named "2D_Visual.ds" and,

- an empty folder named "Clustering_Result".


===================================================================================
 Dataset Format
===================================================================================
The dataset should be given in a text file of the following format:

- In the first line, there are 2 integers: the number $n$ of points in this dataset and the dimensionality $d$ of these points.

- In every subsequent line, there are $d+1$ integers: the id of the point and the $d$ coordinate values, where the id is an integer from 1 to $n$ and each coordinate value is an integer in the range of [0, 10^5].

For instance, the first 6 lines of the sample dataset "2D_Visula.ds" are as shown below:
1000	2
1	10133	92066
2	7133	92600
3	8666	91666
4	7333	90266
5	7933	92000

The first line indicates that there are 1000 2-dimensional points in this dataset. In the second line, there are 3 integers: 1 is the id of this point and (10133, 92066) are its cooridnates. Namely, this line specifies that the point with id = 1 is (10133, 92066). Analogously, the rest four lines above specify the coordinates of the point with id = 2, 3, 4, 5, respectively.  


=================================================================================== 
 Command Line Usage
===================================================================================
Options:
-algo	{integer}	 the algorithm you choose
-n	{integer}	 the number of points
-d	{integer}	 the dimensionality of points
-ds	{string}	 the file path of the dataset
-rf	{string}	 the folder path of the clustering result
-r	{double}	 the radius Eps for DBSCAN
-k	{double}	 the core point threshold MinPts for DBSCAN
-rho	{double}	 the approximation for $rho$-Approximate DBSCAN

Additional Option for our ExactDBSCAN and ApprDBSCAN:
-prg	{0, 1}	 	 If prg = 1, the algorithm runs faster in a progressive manner but consumes more memory space. In default, prg = 0. 

** Note that in our experiment, we always set prg = 0. Namely, this option is always off. **

0. Run Seed Spreader Algorithm to generate random data.
   Parameter list: -algo 0 -n -d -ds

1. Run our exact Algorithm.
   Parameter list: -algo 1 -n -d -r -k -ds -rf [-prg 0]

2. Run our approximate Algorithm.
   Parameter list: -algo 2 -n -d -r -k -ds -rf -rho [-prg 0]

3. Run the GriDBSCAN Algorithm.
   Parameter list: -algo 3 -n -d -r -k -ds -rf

4. Run the Original DBSCAN Algorithm with R-tree.
   Parameter list: -algo 4 -n -d -r -k -ds -rf

5. Run the Original DBSCAN Algorithm with brute-force.
   Parameter list: -algo 5 -n -d -r -k -ds -rf

Let us take the DBSCAN clutering with eps = 5000 and MinPts = 20 on "2D_Visual.ds" as an example. 

Firstly, change the directory to the unzipped folder.

The command line for running our exact DBSCAN algorithm is:

./DBSCAN -algo 1 -n 1000 -d 2 -r 5000 -k 20 -ds "./2D_Visual.ds" -rf "./Clustering_Result/" 

(The command lines for running GriDBSCAN, Original DBSCAN with R-tree and Original DBSCAN with brute-force are the same as above but setting -algo to be 3, 4 and 5, respectively.)

The command line for running our approximate DBSCAN algorithm with rho = 0.001 is :

./DBSCAN -algo 2 -n 1000 -d 2 -r 5000 -k 20 -ds "./2D_Visual.ds" -rf "./Clustering_Result/" -rho 0.001


===================================================================================
 Output Format
===================================================================================
The output of our program consists of two parts:

- In the terminal, the running time of the algorithm is output.

- In the result folder (the folder "./Clustering_Result/" in the previous example), a folder named by the algorithm is generated (if it does not exist) which stores the clustering results. More specifically, each cluster is stored in a single file in text format. In each such file, the first line specifies the number of points in the corresponding cluster, and the other lines specify the points of the cluster in the form of: id and coordinates (the same format as the points  in the dataset).

The name of the folder generated for different algorithms are as follows:

*Folder_Name*		*Algorithm* 
OurExact		Our exact algorithm (-algo 1)
OurAppr			Our approximate algorithm (-algo 2)
GriDBSCAN		GriDBSCAN (-algo 3)
OrgDBSCAN_Rtree		The original DBSCAN algorithm with R-tree (-algo 4)
OrgDBSCAN_BF		The original DBSCAN algorithm with brute-force (-algo 5)

** Please leave these folders *EMPTY* before you run the program. Otherwise, the new clutering results will overwrite (part of) the old ones. **

For example, run the command line

./DBSCAN -algo 1 -n 1000 -d 2 -r 5000 -k 20 -ds "./2D_Visual.ds" -rf "./Clustering_Result/" 

A folder named "OurExact" under the path "./Clustering_Result/" is generated. In the "OurExact" folder, there are 4 files named: "Cluster_1", "Cluster_2", "Cluster_3" and "Cluster_4" meaning that there are 4 clusters.

The first 6 lines in the file (e.g., "Cluster_1") are as shown below:
274
1	10133	92066
7	13066	92466
8	11333	93600
10	9066	94666
11	13466	95400

The number 274 in the first line means that there are 274 points in this cluster. The second line means that the point with id = 1 whose coordinates are (10133, 92066) is in this cluster. Anologously, the rest lines specify which points are in this cluster.

  






