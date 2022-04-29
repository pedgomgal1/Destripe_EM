# -*- coding: utf-8 -*-
# """
# Created on Thu Apr 21 15:22:53 2022

# @author: Pedro Gomez Galvez
# """
#matplotlib inline
from matplotlib import pyplot as plt

from skimage.io import imread                                                   
import numpy as np
from myFunctions_destripe import main_mask_fft
#import cv2


imagePath=r"C:\Users\Pedro Gomez Galvez\Documents\Lab\Destripe_EM\FIBdeSEMAna_cropped_8bits_resize.tif"

img_noisy=imread(imagePath);#.astype(np.double);
#img_noisy = cv2.resize(img_noisy, dsize=(1064,898), interpolation=cv2.INTER_LINEAR)


plt.imshow(img_noisy, cmap='gray')


print("-----------------------------------------------\n")
print("Dataset Loaded\n")
print("-----------------------------------------------\n")

shpImg=img_noisy.shape
mask_fft=np.zeros(shpImg, dtype=float);

print('-----------------------------------------------\n');
print('Starting Mask detection\n');
print('-----------------------------------------------\n');

mask_fft = main_mask_fft(img_noisy);
    
# print('-----------------------------------------------\n');
# print('Mask Detection Done\n');
# print("-----------------------------------------------\n")

# %%
# fprintf('-----------------------------------------------\n');
# fprintf('Simple approach recovery using Dykstra\n');
# fprintf('-----------------------------------------------\n');
# Rec_Zero_fft = zeros(m,n,k);
# for index = 1:k
# 	Rec_Zero_fft(:,:,index) = 255*FFT_Zero_proj(img_noisy(:,:,index),mask_fft(:,:,index));
# 	fprintf('Slice %4d of total %4d is done\n',index,k);
# end



    
