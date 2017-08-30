# SliVo

Slice to volume registration pipeline developed for 2D histological slide to 3D micro computed tomography registration. Details of the algorithm are explained in 

Natalia Chicherova, Ketut Fundana, Bert Müller, Philippe C. Cattin,
[Histology to μCT Data Matching Using Landmarks and a Density Biased RANSAC](https://link.springer.com/chapter/10.1007/978-3-319-10404-1_31), Lecture Notes in Computer Science - MICCAI 2014 8673: 243–250.

# Example
```Matlab
% example of use to register a 2D slice to a 3D volume
clear all
close all
% load 3D volume 
volume_dir = './Data/';
filename=[volume_dir 'CT_data.mat'];
load(filename);

% load grayscale histology
histology_dir = './Data/';
filename_hist = [histology_dir 'Histology.mat'];
load(filename_hist);

% localize a histological slide in a 3D data
registerSliceToVolume(CT_data,Histology,'lower_limit',100, 'upper_limit', 890,'calculate_features', 1);
```

# Implementation details
- histology must be converted to grayscale ``rgb2gray(image)``
- better to use 8bit images
- better limit or cut the 3D volume images with low ROI ``'upper_limit','lower_limit'``.  


# Acknowledgements
1. D. Kroon for [SURF](http://ch.mathworks.com/matlabcentral/fileexchange/28300-opensurf--including-image-warp-)
2. Peter Kovesi for [RANSAC](http://www.peterkovesi.com/matlabfns/)