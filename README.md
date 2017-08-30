# SliVo

Slice to volume registration pipeline developed for 2D histological slide to 3D micro computed tomography registration. Details of the algorithm are explained in 

Natalia Chicherova, Ketut Fundana, Bert Müller, Philippe C. Cattin,
[Histology to μCT Data Matching Using Landmarks and a Density Biased RANSAC](https://link.springer.com/chapter/10.1007/978-3-319-10404-1_31), Lecture Notes in Computer Science - MICCAI 2014 8673: 243–250.

# Implementation details
- histology must be converted to grayscale ``rgb2gray(image)``
- better to use 8bit images
- better limit or cut the 3D volume images with low ROI ``'upper_limit','lower_limit'``.  


# Acknowledgements
1. D. Kroon for SURF
2. Peter Kovesi for RANSAC