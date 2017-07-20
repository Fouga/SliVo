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
% Copyright © 2017 Natalia Chicherova.

if ~isempty(v)
    if iscell(v{1})
        v = v{1};
    end
end

%% options that may be specified by the user
options.feature_detector           = []; % default SURF
options.filter_radius_ratio             = 2.8;   % default value = 2.8
options.lower_limit             = [];   % default value = 
options.upper_limit    = [];   % default value = 
options.angle             = pi/8;   % default value = pi/8
options.calculate_features      = 1; % if you do not need to calculate features
options.number             = 0;   % number of points to leave after filtering, by default it is calculated automatically

options.folder_source           = [];
options.folder_destination      = [];
options.folder_features         = [];
options.folder_matches          = [];
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
    else
        warning('SliVo:parseInput', 'Expected numeric value\n');
    end
end


function param = getFolder(v,i)

param = [];
if i+1 <= numel(v)
    if isdir(v{i+1})
        param = v{i+1};
    else
        if isstr(v{i+1});
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


