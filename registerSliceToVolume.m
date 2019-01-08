function found_slice = registerSliceToVolume(Data3D,Slice2D,varargin)
% Registation pipeline for localizing an input 2D image (Slice2D) in a 3D
% volume (Data3D). Includes 3 main steps: find corresponding features between
% Slice2D and each image in the Data3D, fit a plane to the highest density 
% of the feature cloud and, finally, optimizing plane location in 3D.  
% 
% Usage:  found_slice = registerSliceToVolume(Data3D,Slice2D,varargin)
%
% Input: Data3D - is a XxYxZ 3D uint8 array. The image volume data should be
%                 passed as a stack of 2D images 
%       Slice2D - is a 2D uint8 image which you want to localize in the
%                 Data3D, for example, histological slice.
%      varargin - (optional) is a structure containing various parameter values needed by
%                slice to volume registration. For more details about the
%                pissible registration parameters, see - SliVo_parseInputs.
%
% Output: found_slice - is a 2D image extracted from the Data3D which plane parameters 
%                       correspond to the best localization of the input Slice2D in 3D,
%                       i.e. it is a position of the input Slice2D in the 3D data (Data3D). 
%         The output results are saved in the otput folder: options.folder_destination
%
% See also: find_features, fit_plane, optimize_plane
%
% From the project, (https://github.com/Fouga/).
% Author: 2019 Natalia Chicherova.
%
% If you use this code please cite my paper:
% Chicherova, N., Hieber, S.E., Khimchenko, A., Bikis, C., Mü ller, B.,
% Cattin, P.C., 2018. Automatic deformable registration of histological
% slides to μCT volume data. Journal of Microscopy 271, 49–61.

if nargin==0
	help(mfilename)
	return
end

% add necessary paths
if (~isdeployed)
    p = mfilename('fullpath');
    p = p(1:end-22);
    addpath(genpath(fullfile(p, 'main'))); 
    addpath(genpath(fullfile(p, 'Feature_detection'))); 
    addpath(genpath(fullfile(p, 'MatlabFns')));
end

options = SliVo_parseInputs(varargin);

options = find_features(Data3D,Slice2D,options);
% display

% fit a plane to the highest feature density
found_slice =fit_plane(Data3D,Slice2D,options);

if options.optimization == true
   found_slice = optimize_plane(Data3D,Slice2D,found_slice,options);
end