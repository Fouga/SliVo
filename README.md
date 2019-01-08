# SliVo

This 2D to 3D registration pipeline is developed to help localizing a histological slice in a micro computed tomography (CT) volume. The algorithms is based on coarse-to-fine approach. First, the position of the 2D histological slice is initialized as a plane in the 3D CT volume and, then, the plane coordinates are optimized using more sensitive similarity measure. 

The extension of the algorithm can be also applied to MRI. The best performance in this case showed not invariant to rotation extended [SS](https://ieeexplore.ieee.org/abstract/document/4270223) descriptor that is called SL1. 

<img src="![pipeline](https://user-images.githubusercontent.com/17926378/50838463-55de5b00-135e-11e9-808a-2f8d29446d97.png)" />


Details of the algorithm are explained in the following papers:

[1] Natalia Chicherova, Ketut Fundana, Bert Müller, Philippe C. Cattin,
[Histology to μCT Data Matching Using Landmarks and a Density Biased RANSAC](https://link.springer.com/chapter/10.1007/978-3-319-10404-1_31), Lecture Notes in Computer Science - MICCAI 2014 8673: 243–250.

[2] Chicherova, N., Hieber, S.E., Khimchenko, A., Bikis, C., Müller, B., Cattin, P.C.
[Automatic deformable registration of histological slides to μCT volume data](https://onlinelibrary.wiley.com/doi/full/10.1111/jmi.12692), Journal of Microscopy 271, 49–61, 2018.

or see in ./pdf/Histology_to_mCT_Data_Matching_using_Landmarks_and_a_Density_Biased_RANSAC.pdf

# Example
```Matlab
% example of use to register a 2D slice to a 3D volume
clear all
close all
% load 3D CT volume 
volume_dir = './Data/3dCT/';
imds =imageDatastore(fullfile(volume_dir), 'FileExtensions', {'.tif'});
[Xsize Ysize] = size(readimage(imds,1));
Zsize = size(imds.Files,1);
CTvolume = zeros(Xsize, Ysize, Zsize, 'uint8');
parfor i = 1:Zsize
    CTvolume(:,:,i) = readimage(imds,i);
end

% load grayscale histology
histologyName = './Data/Resized_histology.png';
Histology = imread(histologyName);

% localize a histological slide in a 3D data
registerSliceToVolume(CTvolume,Histology,'lower_limit',40, 'upper_limit', 400,...
    'calculate_features', 1, 'optimization',1,'angle',pi/18,'rotation_invariance',1,...
    'method2dregistration','NCC_NMI', 'radius',2.2);
```

# Implementation details
- histology must be converted to grayscale ``rgb2gray(image)``
- better to use 8bit images
- better limit or cut the 3D volume images with low ROI ``'upper_limit','lower_limit'``.  
- the histology and the 2D slice obtained after RANSAC fit should be as precise as possible registered in 2D. One can use various 2D-2D registration algorithms to improve the performance. The correspondences based Ransac Homography registration (as in [1]((https://link.springer.com/chapter/10.1007/978-3-319-10404-1_31))) is not implemented. 


# Acknowledgements
1. D. Kroon for [SURF](http://ch.mathworks.com/matlabcentral/fileexchange/28300-opensurf--including-image-warp-)
2. Peter Kovesi for [RANSAC](http://www.peterkovesi.com/matlabfns/)
