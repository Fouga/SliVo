% example of use to register a 2D slice to a 3D volume
clear all
close all
% load 3D volume 
volume_dir = '/media/natasha/LMP_DATA/Natasha/Histology_registration_final/Data/';
filename=[volume_dir '13_301x301x969_cut.mat'];
load(filename);

% load histology
histology_dir = '/media/natasha/LMP_DATA/Natasha/Histology_registration_final/Data/';
filename_hist = [histology_dir 'Cut_Histology4.mat'];
load(filename_hist);

% tmp
registerSliceToVolume(data_cyl_cut,Hist,'lower_limit',100, 'upper_limit', 890,'calculate_features', 1);



