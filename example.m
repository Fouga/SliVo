% example of use to register a 2D slice to a 3D volume

% load 3D volume 
volume_dir = './Data/';
filename=[volume_dir '13_301x301x969_cut.mat'];
load(filename);

% load histology
histology_dir = './Data/';
filename_hist = [histology_dir 'Cut_Histology1.mat'];
load(filename_hist);

% tmp
registerSliceToVolume(data_cyl_cut,Hist,'lower_limit',100, 'upper_limit', 890,'calculate_features', 1);

% registerSliceToVolume(Data3D,histology);

