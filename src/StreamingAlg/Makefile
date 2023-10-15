biodeg:
	python ./main.py \
		--data=../../datasets/biodeg_data.csv \
		--buffer_size=50 \
		--rho=0.5\
		--r=5 \
		--minpts=10 \
		--label_out=./output/biodeg_label.csv
	python ./Analysis.py \
		--y=../../datasets/biodeg_label.csv \
		--label=./output/biodeg_label.csv

mnist:
	python ./DBSCAN.py \
		--data=../../datasets/mnist_aug_small_data.csv \
		--buffer_size=100 \
		--rho=0.5 \
		--eps=200.0 \
		--minpts=10 \
		--label_out=./output/mnist_aug_small_label.csv
	python ./Analysis.py \
		--y=../../datasets/mnist_aug_small_label \
		--label=./output/mnist_aug_small_label.csv

cifar:
	python ./DBSCAN.py \
		--data=../../datasets/cifar_aug_small_data.csv \
		--buffer_size=100 \
		--rho=1 \
		--eps=5000.0 \
		--minpts=10 \
		--label_out=./output/cifar_aug_small_label.csv
	python ./Analysis.py \
		--y=../../datasets/cifar_aug_small_label \
		--label=./output/cifar_aug_small_label.csv

fashion_mnist:
	python ./DBSCAN.py \
		--data=../../datasets/fashion_mnist_aug_small_data.csv \
		--buffer_size=100 \
		--rho=0.5 \
		--eps=8000.0 \
		--minpts=10 \
		--label_out=./output/fashion_mnist_aug_small_label.csv
	python ./Analysis.py \
		--y=../../datasets/fashion_mnist_aug_small_label \
		--label=./output/fashion_mnist_aug_small_label.csv

uspshw:
	python ./DBSCAN.py \
		--data=../../datasets/uspshw_aug_small_data.csv \
		--buffer_size=100 \
		--rho=0.5 \
		--eps=500.0 \
		--minpts=10 \
		--label_out=./output/uspshw_aug_small_label.csv
	python ./Analysis.py \
		--y=../../datasets/uspshw_aug_small_label \
		--label=./output/uspshw_aug_small_label.csv

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

cancer:
	python ./main.py \
		--data=../../datasets/breastcancer_data.csv \
		--buffer_size=100 \
		--rho=0.5 \
		--r=30 \
		--minpts=4 \
		--label_out=./output/breastcancer_label.csv
	python ./Analysis.py \
		--y=../../datasets/breastcancer_label.csv \
		--label=./output/breastcancer_label.csv

arrhythmia:
	python ./main.py \
		--data=../../datasets/arrhythmia_data.csv \
		--buffer_size=100 \
		--rho=0.5\
		--r=200 \
		--minpts=10 \
		--label_out=./output/arrhythmia_label.csv
	python ./Analysis.py \
		--y=../../datasets/arrhythmia_label.csv \
		--label=./output/arrhythmia_label.csv