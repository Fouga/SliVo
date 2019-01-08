function Slice2D_optimal = optimize_plane(Data3D,Slice2D,found_slice,options)
% Optimizes the localization of the found plane with Ransac. The found with fit_plane 
% plane parameters serve as initialization for the optimization algortihm for non-smooth 
% cost function (NMI). For simplicity, the normal vector to the found plane 
% is converted to spherical coordinate system. 
% The space of initial plane paramaters is sampled 20 times around the one
% found with Ransac and finally the plane with the highest NMI is chosen.
%
% In order to use NMI for plane paramaters optimization the 2D slice should
% be very well aligned in 2D with the 3D data. The rotation, scaling and
% shift should be found either automatically or manually (if the input is already aligned). 
%
% Usage:  Slice2D_optimal = optimize_plane(Data3D,Slice2D,found_slice,options)
%
% Input: Data3D - is a XxYxZ 3D uint8 array. The image volume data should be
%                 passed as a stack of 2D images 
%       Slice2D - is a 2D uint8 image which you want to localize in the
%                 Data3D, for example, histological slice.
%   found_slice - is a 2D image from the Data3D found as a fit to the highest density of
%                 3D correspondences cloud. 
%       options - is a structure containing various parameter values needed by
%                slice to volume registration 
%
% Output: Slice2D_optimal - is a 2D image extracted from the Data3D which plane parameters 
%                           correspond to the maximum of NMI between Slice2D and Data3D. 
%         The output results are saved in the otput folder: options.folder_destination
%
% See also: fit_plane, register2D2D, fminsearch
%
%
% From the project, (https://github.com/Fouga/).
% Author: 2019 Natalia Chicherova.
%
% If you use this code please cite my paper:
% Chicherova, N., Hieber, S.E., Khimchenko, A., Bikis, C., Mueler, B.,
% Cattin, P.C., 2018. Automatic deformable registration of histological
% slides to μCT volume data. Journal of Microscopy 271, 49–61.

Slice2D_registered2D = register2D2D(Slice2D,found_slice,options);

Slice2D_optimal = optimize_spherical_coor(Data3D,Slice2D_registered2D,options);






function Slice2D_optimal = optimize_spherical_coor(Data3D,Slice2D_registered2D,options)


if ~exist(fullfile(options.folder_destination, 'Optimized_slice'))
    mkdir(fullfile(options.folder_destination, 'Optimized_slice'))
end
options.folder_optimal = fullfile(options.folder_destination, 'Optimized_slice');  
% set angle limits to the optimization in radians
max_angle =options.angle;

[X_size Y_size Z_size] = size(Data3D);
[x_mesh, y_mesh]=meshgrid(1:Y_size,1:X_size);


%% Multisearch fminsearch
norm_al = options.angle; norm_phi = 2*pi; norm_d = Z_size;
num = 20; % number of random trials for optimization
Optim_randSearch = zeros(num+2,4);
x2_rnd = (0 + (360-0)*rand(num,1))*pi/180/norm_phi;
x1_rnd = (0 + (max_angle-0)*rand(num,1))*pi/180/norm_al;

initial = 'Plane parameters optimization';
h = waitbar(0,initial);
m = 1;
myfilter = fspecial('gaussian',[3 3], 0.5);
hist_smooth = imfilter(Slice2D_registered2D, myfilter, 'replicate','same');
    
% load normal vector
filename_B = fullfile(options.folder_matches, 'Normal_vec_toPlane.mat');
load(filename_B,'B');
% convert to spherical coordinates
    l = B(1);
    B(1) = B(2);
    B(2) = l;
    B_unit = B./sqrt(B(1)^2 + B(2)^2 + B(3)^2);
    [ alpha, phi, d ] = planeParams( B_unit ); % in degree
    % change to cartesian initial guess
    if phi < 0
        phi = 2*pi + phi;
    end
    sph = [ alpha, phi, d ];
    save (fullfile(options.folder_matches, 'Normal_vec_toPlane_polar.mat'),'sph','-mat');
    
    
x0 = [sph(1)/norm_al, sph(2)/norm_phi, sph(3)/norm_d]';%normalise
d_lb = (sph(3)-80)/norm_d;% the search in Z is limited by 80 slices up and down
d_ub = (sph(3)+80)/norm_d;
    
% put a NMI value of the slice found with Ransac with no optimization step
fval = optim_func(x0,Data3D, hist_smooth,max_angle);
Optim_randSearch(m,:) = [sph, fval] ;
m = m+1;
        
% Optimisation from initial RANSAC guess
LB = [-max_angle/norm_al -2*pi/norm_phi d_lb];% normalized optimization constraints 
UB = [max_angle/norm_al 2*pi/norm_phi d_ub];
func_handle = @(x0)optim_func(x0,Data3D,hist_smooth,max_angle);
[x0_search, fval] = fminsearchbnd(func_handle, x0,LB,UB);

x0_opt = [x0_search(1)*norm_al, x0_search(2)*norm_phi, x0_search(3)*norm_d];
Optim_randSearch(m,:) = [x0_opt, fval];
m= m+1;
        
%%%% Multi random start 
x3_rnd = (d_lb + (d_ub - d_lb)*rand(num,1));

for i = 1:num
    mess = sprintf('%i itteration of the plane parameters optimization', i);
    waitbar(i/num,h,mess);
    %             try_values = [];
    x0 = [x1_rnd(i) x2_rnd(i) x3_rnd(i)];
    % convert to cartesian to find plane
    alpha = x0(1)*norm_al; phi = x0(2)*norm_phi; d = x0(3)*norm_d;
    len = 1.000;
    B_ = len * sin(alpha)*sin(phi);
    C = len * cos(alpha);
    A = len * sin(alpha)*cos(phi);
    D = -d * C;
    Z_LB =min(min(-(y_mesh.*B_+ x_mesh.*A + D)/C));
    Z_UB = max(max(-(y_mesh.*B_+ x_mesh.*A + D)/C));
    w = 1;
    while (Z_LB < 1  && w <1000) || (Z_UB > Z_size && w < 1000)
        x2_new = (0 + (360-0)*rand(1,1))*pi/180/norm_phi;
        x1_new = (0 + (alpha_max-0)*rand(1,1))*pi/180/norm_al;
        x3_new = (d_lb + (d_ub - d_lb)*rand(1,1));
        x0 = [x1_new x2_new x3_new];
        alpha = x0(1)*norm_al; phi = x0(2)*norm_phi; d = x0(3)*norm_d;
        len = 1.000;
        B_ = len * sin(alpha)*sin(phi);
        C = len * cos(alpha);
        A = len * sin(alpha)*cos(phi);
        D = -d * C;
        Z_LB =min(min(-(y_mesh.*B_+ x_mesh.*A + D)/C));
        Z_UB = max(max(-(y_mesh.*B_+ x_mesh.*A + D)/C));
        w = w+1;
        if w == 999
            display('did not find good random initialization!');
            break;
        end
    end
            
    func_handle = @(x0)optim_func(x0,Data3D,hist_smooth,max_angle);
    [x0_search,fval] = fminsearchbnd(func_handle, x0,LB,UB);
    x0_opt = [x0_search(1)*norm_al, x0_search(2)*norm_phi, x0_search(3)*norm_d];
    Optim_randSearch(m,:) = [x0_opt, fval];
    m = m+1;
end
close(h);

save(fullfile(options.folder_optimal, 'OptimalRandSearchCoor.mat'),'Optim_randSearch'); 

% find plane parameters with the lowest cost

[val_min ind_min] = min(Optim_randSearch(:,end));% find min NMI
Sph_optimal = Optim_randSearch(ind_min,1:3);

save(fullfile(options.folder_optimal, 'NormalVec_optimal_polar.mat'),'Sph_optimal'); 

%% plot results

% optimisation slice
x0_min = Sph_optimal;
alpha = x0_min(1); phi = x0_min(2); d = x0_min(3);
len = 1.000;
B_ = len * sin(alpha)*sin(phi);
C = len * cos(alpha);
A = len * sin(alpha)*cos(phi);
D = -d * C;
Z_opt =-(y_mesh.*B_+ x_mesh.*A + D)/C;
Slice2D_optimal = interp3(single(Data3D),x_mesh,y_mesh,Z_opt);

% ransac slice
ransac_coor = Optim_randSearch(1,:);
alpha = ransac_coor(1); phi = ransac_coor(2); d = ransac_coor(3);
len = 1.000;
B_ = len * sin(alpha)*sin(phi);
C = len * cos(alpha);
A = len * sin(alpha)*cos(phi);
D = -d * C;
Z_Ransac =-(y_mesh.*B_+ x_mesh.*A + D)/C;  
Slice_Ransac = interp3(single(Data3D),x_mesh,y_mesh,Z_Ransac); 

close all
figure, subplot(2,3,1), imshow(Slice2D_optimal,[]), title('Found optim slice')
        subplot(2,3,2), imshow(Slice2D_registered2D,[]), title('Input slice')
        subplot(2,3,3), imshow(single(Slice2D_registered2D) - Slice2D_optimal,[]), title('Difference Found optim slice')

        subplot(2,3,4), imshow(Slice_Ransac,[]), title('Ransac')
        subplot(2,3,5), imshow(Slice2D_registered2D,[]), title('Input slice')
        subplot(2,3,6), imshow(single(Slice2D_registered2D) - Slice_Ransac,[]), title('Difference Ransac')

print ('-dtiff', fullfile(options.folder_destination, 'Slice_comparison_image.tif'));
imwrite(uint8(Slice2D_optimal), fullfile(options.folder_destination, 'Final_registration_result.tif'));


  