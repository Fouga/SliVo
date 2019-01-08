function options = SliVo_parseInputs(v)  
% Parses the arguments and fills the options structure with
% default parameter values 
%
% Usage:  OPTIONS = SliVo_parseInputs(v)
%
% Input:  V a cell containing arguments from the registerSliceToVolume (from varargin)
%
% Output: OPTIONS a structure containing various parameter values needed by
%         slice to volume registration
%
% See also: registerSliceToVolume
%
% From the project, (https://github.com/Fouga/).
% Author: 2019 Natalia Chicherova.
%
% If you use this code please cite my paper:
% Chicherova, N., Hieber, S.E., Khimchenko, A., Bikis, C., Müuller, B.,
% Cattin, P.C., 2018. Automatic deformable registration of histological
% slides to μCT volume data. Journal of Microscopy 271, 49–61.


if ~isempty(v)
    if iscell(v{1})
        v = v{1};
    end
end

%% options that may be specified by the user
options.feature_detector           = 'SURF'; % default SURF
options.rotation_invariance           = 1; % use rotation invariant Feature detectro by default
options.filter_radius_ratio             = 2.8;   % default value = 2.8
options.lower_limit             = [];   % default value = 
options.upper_limit    = [];   % default value = 
options.angle             = pi/10;   % default value = pi/10
options.calculate_features      = 1; % if you do not need to calculate features
options.number             = 0;   % number of points to leave after filtering, by default it is calculated automatically
options.optimization       = 1; % after the initilization with the feature detectors,
% the plane parameters can be optimized using NMI. By default the
% optimization is done
options.method2dregistration = 'NCC_NMI'; % method for 2D 2D registration

options.folder_source           = [];
options.folder_destination      = [];
options.folder_features         = [];
options.folder_matches          = [];
options.folder_optimal          = [];
options.filenames               = {};

if nargin==0 || isempty(v)
    return
end

% handle the input parameters, provided in (string, param) pairs
for i = 1:numel(v)
    if ischar(v{i})
            switch lower(v{i})
            case {'filter_radius_ratio', 'radius'}
                options.filter_radius_ratio = getParam(v,i);            
             case 'lower_limit'
                options.lower_limit   = getParam(v,i);   
             case 'upper_limit'
                options.upper_limit   = getParam(v,i);         
             case 'calculate_features'
                options.calculate_features   = getParam(v,i);     
             case 'angle'
                options.angle   = getParam(v,i);  
             case 'number'
                options.number   = getParam(v,i);  
             case 'feature_detector'
                options.feature_detector  = getParam(v,i);  
             case 'rotation_invariance'                
                options.rotation_invariance = getParam(v,i);  
             case 'optimization'
                options.optimization  = getParam(v,i);  
             case 'method2dregistration'
                options.method2dregistration  = getParam(v,i);  
             case 'source'
                options.folder_source  = getFolder(v,i);  
            case 'destination'
                options.folder_destination = getFolder(v,i);   
            end
    end 
end

    % sort the options alphabetically so they are easier to read
options = orderfields(options);



function param = getParam(v,i)

param = [];
if i+1 <= numel(v)
    if isnumeric(v{i+1})
        param = v{i+1};
    elseif ischar(v{i+1})
           param = v{i+1};
    else
        warning('SliVo:parseInput', 'The parameter is not correct\n');
    end
end


function param = getFolder(v,i)

param = [];
if i+1 <= numel(v)
    if isfolder(v{i+1})
        param = v{i+1};
    else
        if ischar(v{i+1})
            param = v{i+1};
            [success message] = mkdir(param);
            if ~success
                warning('SliVo:parseInput', message);
            end
        else
            warning('SliVo:parseInput', 'Expected destination path value\n');
        end
    end
end


