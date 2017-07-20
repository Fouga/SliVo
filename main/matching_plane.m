function [B,inliers]=matching_plane(surf_3D, options)
%%%                             THE CODE

X_size = options.size(1);
Y_size = options.size(2);
radius = X_size/options.filter_radius_ratio;

%Build 3D data set out of coordinates
data_3D = [];
for i=options.lower_limit:options.upper_limit 
    n = size(surf_3D{i},1);
    data_3D = [data_3D; surf_3D{i} repmat(i,n,1)];   
end

matr = zeros(X_size,Y_size,options.upper_limit ,'uint8');%(row, column, )

%put ones where it has a feature
l = 1;
for j=options.lower_limit:options.upper_limit 
    sliceData = data_3D(:,3) == j; %take all Z rows which =j
    sliceX = data_3D(sliceData,2);%column
    sliceY = data_3D(sliceData,1);
    dlina = length(sliceX);
    Z = j;
    for i = 1:dlina
         X = round(sliceX(i));
         Y = round(sliceY(i));
         if ((X-X_size/2)^2 +(Y-Y_size/2)^2) < radius^2 
             matr(X,Y,Z) = matr(X,Y,Z) + 1;
             Rounded_coor(l,1:3) = [Y X Z];
             l = l + 1;
         end
    end      
end

%first smooth the data to find locations with high densities
siz=181;%181; %kernal size
sigma = siz/2.34/2;%smoothing size
siz = (siz-1)./2;
x = linspace(-siz, siz, siz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
matr_smoothX = convn(matr,gaussFilter,'same');
clear matr
matr_smoothXY = convn(matr_smoothX,gaussFilter','same');
clear matr_smoothX
siz=5;
sigma = siz/2.34/2;%sigma = 0.65;
siz = (siz-1)./2;
x = linspace(-siz, siz, siz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); 
z=zeros(1,1,numel(x));
z(:,:,:)=gaussFilter; 
matr_smoothXYZ = convn(matr_smoothXY,z,'same');
clear matr_smoothXY

%save new built coordinates
%Rounded_coor=round(data_3D);
X_r=Rounded_coor(:,2);
Y_r=Rounded_coor(:,1);
Z_r=Rounded_coor(:,3);
clear Rounded_coor
%read the values in coordiantes where we put one
for k=1:length(X_r)
    pointX=X_r(k);
    pointY=Y_r(k);
    pointZ=Z_r(k);
    v(k)=matr_smoothXYZ(pointX,pointY,pointZ);
end
clear matr_smoothXYZ
%sort them and pick first 2000
[vsort,idx] = sort(v,'ascend');%idx-where a big value stands in v
if options.number ==0
    total_num = length(X_r)
    num = find_equi_number(total_num)
else num = options.number;
end

% if length(X_r) > num
    num_retain = num; 
    x_ransac = X_r(idx(end-num_retain:end));
    y_ransac = Y_r(idx(end-num_retain:end));
    z_ransac = Z_r(idx(end-num_retain:end));
    plane_ransac=[x_ransac,y_ransac,z_ransac];
    weights=vsort(end-num_retain:end);
% else 
%     x_ransac = X_r;
%     y_ransac = Y_r;
%     z_ransac = Z_r;
%     plane_ransac=[x_ransac,y_ransac,z_ransac];
%     weights=vsort;
% end

%%%%%%%%%%%%%%%%%ransac fit a plane to the picked points
t=10;
figure,
x = plane_ransac(:, 1);
y = plane_ransac(:, 2);
z = plane_ransac(:, 3);
plot3(x,y,z,'*b');

% weights=vsort(end-num_retain:end);
    [B, P, inliers] = ransacfitplane_density(plane_ransac',t,0,weights,options.angle);

hold on;
[xmesh, ymesh] = meshgrid(linspace(min(x), max(x), 10), linspace(min(y), max(y), 10));
z=-(ymesh.*B(2) + xmesh.*B(1) + B(4))/B(3);
surf(xmesh, ymesh, z, 'facecolor', 'none');
hold off

% figure,
% x = X_r;
% y = Y_r;
% z = Z_r;
% plot3(x,y,z,'*b');

