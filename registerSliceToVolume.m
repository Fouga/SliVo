function found_slice = registerSliceToVolume(Data3D,Slice2D,varargin)

if nargin==0
	help(mfilename)
	return
end

% add necessary paths
if (~isdeployed)
    p = mfilename('fullpath');
    p = p(1:end-22);
    addpath([p '/main/']); 
    addpath(genpath([p '/Feature_detection/'])); 
    addpath(genpath([p '/MatlabFns/']));
end

options = SliVo_parseInputs(varargin);

options = find_features(Data3D,Slice2D,options);

% fit a plane to the highest feature density
found_slice =fit_plane(Data3D,Slice2D,options);

