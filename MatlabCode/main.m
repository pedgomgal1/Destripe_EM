%% start
clc; clear; close all;
addpath(genpath(pwd));
%% Load dataset 3D

%load(fullfile('data','dataset_3520_6580_5600_rot.mat'));
img = double(imread('../../FIBdeSEMAna_cropped_8bits.tif'));
% img=imresize(img,0.1);
slices=img;
%slices = double(slices(1:end,1:end,1:2:20));
img_noisy = slices/255;
fprintf('-----------------------------------------------\n');
fprintf('Dataset Loaded\n');
fprintf('-----------------------------------------------\n');

[m,n,k] = size(slices);
mask_fft = zeros(m,n,k);
fprintf('-----------------------------------------------\n');
fprintf('Starting Mask detection\n');
fprintf('-----------------------------------------------\n');
for index = 1:k
	img = squeeze(slices(:,:,index));
	mask_fft(:,:,index) = main_mask_fft(img);
	fprintf('Slice %4d of total %4d is done\n',index,k);
end

fprintf('-----------------------------------------------\n');
fprintf('Mask Detection Done\n');
fprintf('-----------------------------------------------\n');

%%
fprintf('-----------------------------------------------\n');
fprintf('Simple approach recovery using Dykstra\n');
fprintf('-----------------------------------------------\n');
Rec_Zero_fft = zeros(m,n,k);
for index = 1:k
	Rec_Zero_fft(:,:,index) = 255*FFT_Zero_proj(img_noisy(:,:,index),mask_fft(:,:,index));
	fprintf('Slice %4d of total %4d is done\n',index,k);
end
%%
fprintf('-----------------------------------------------\n');
fprintf('Complecated approach using 2D FFT and 3D TV\n');
fprintf('-----------------------------------------------\n');

Rec_FFT_TV = 255*FFT2_TV3_L_Proj(img_noisy,mask_fft);


fprintf('-----------------------------------------------\n');
fprintf('recovery Done!\n');
fprintf('-----------------------------------------------\n');

