minimal_scene
=============

Code for the CVPR 2014 paper "Minimal Scene Descriptions from Structure from Motion Models"

Platform tested: Ubuntu 12.04, Ubuntu 14.04
Packages required: g++, libtbb-dev, zlib1g-dev


Usage: point_reduce <img_pt_mat file> <list keys file> <num images> <num points> <k cover> <percentage> <point idx file> [record dist] [use k-cover] [reduce memory] [ip weight] [min prob] [threshold as prob] [none 0; read mean 1; write mean 2] [mean file] [use binary] [load pt idx] [cdf file] [cdf as prob]

<> means required arguments and [] means optional arguments. Here is a subset of arguments explained:

<img_pt_mat file>: the file that contains the bi-partite graph between points and images. Each line should be in the format of 
        pt_idx1 img_idx1 pos1
        pt_idx2 img_idx2 pos2
        ...
    , where pt_idx is the point index, img_idx is the image index, and pos is the index of which the point's descriptor appear in the feature (e.g. SIFT) list of that image.

<list keys file>: each line of this file should be the location of a SIFT key file. The order should be consistent with image index (0 based). 

<num images>: total number of images

<num points>: total number of points

<k cover>: K used in the probabilistic K-cover algorithm

<percentage>: The target percentage of coverage over all images

<point idx file>: The output file of selected point indices

[record dist]: Record distribution of each image being covered, for debugging purposes. Normally could be set to 0.

[use k-cover]: 1 for running the baseline K-cover algorithm, otherwise 0. 

[reduce memory]: whether try to reduce memory usage, 0 for not to, which is faster.

[ip weight]: The "p" in the paper used as the constant probability of an visibility event. (0.6 used in the paper)

[threshold as prob]: Minimum probability for "covering" an image. (p_{min} in the paper)

Optionally, to save time of reading all the SIFT key files and compute the mean descriptors of all the points, one can choose to run the program once and write the mean descriptors of each point as a file, and all subsequent running instances could choose to read the mean file. [use binary] dictates whether the binary format is used for these operations.

Example command:
point_reduce data/img_pt_pos_mat data/list.key.txt 6044 1886884 12 0.99 selected_idx 0 0 0 0.6 0.99 0 1 006/data/mean_desc.bin 1