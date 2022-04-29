# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:35 2022

@author: Pedro Gomez Galvez
"""

# myFunctions_destripe.py
import math
from collections import namedtuple
import numpy as np
from scipy.stats import rayleigh 
from scipy.signal import medfilt2d
import scipy



def main_mask_fft(img):
    
    options = namedtuple('options', ['allout', 'angle_threshold','random_on','cutoff_threshold','dilation_mode','dilation_size','mask_mode','medfilt_size','box_mask']);
    options.allout=0;
    options.angle_threshold=0.1;
    options.random_on=1;
    options.cutoff_threshold=1e-3;
    options.dilation_mode='square';
    options.dilation_size=5;
    options.mask_mode='medfilt';
    options.medfilt_size=[5,5];
    options.box_mask=False;

    [mask_fft,p_value_mask, statistics, positions, indecies, ring] = get_mask_fft(img,options);
  #   [angle,binom_val] = get_stripe_angle(p_value_mask, ring);
  #   [min_ind,min_val,binom_val_box] = get_stripe_best_box(p_value_mask,ring);
  #   d = min_ind(1);
  #   mask_fft_modif = numpy.fft.ifftshift(mask_fft);
  #   mask_fft_modif(1:ring.center(1)-d-10,:)=0;
  #   mask_fft_modif(ring.center(1)+d+10:end,:)=0;

    
 	# SE = strel('square',options.dialation_size);
 	# mask_fft_modif =  imdilate(mask_fft_modif,SE);

 	
 	
 	# center = math.floor(size(img)/2)+1;
 	# center_box_size = 4;
 	# mask_fft_modif(center(1)-center_box_size:center(1)+center_box_size,center(2)-center_box_size:center(2)+center_box_size)=0;
 	
 	# mask_fft_modif = fftshift(mask_fft_modif);
 	# mask_fft = mask_fft_modif;
    
    return mask_fft

def get_mask_fft (img,options):
    # fft calculation
    im_fft = np.fft.fft2(img);
    img_fft_shift = np.fft.fftshift(im_fft);

    # Ring parameters
    ring = ring_initialize(img_fft_shift.shape);

    # main body
    [p_value_mask] = generate_p_value_mask(img_fft_shift,ring,options);
    binary_mask_medfilt = generate_binary_mask(p_value_mask,options);
    mask_fft = np.fft.ifftshift(binary_mask_medfilt);
    
    return mask_fft,p_value_mask, ring

def ring_initialize (img_size):
    ring=namedtuple('ring',['center','radius_init','halfwidth','radius_end','radius','number'])
    ring.center = [math.floor(img_size[0]/2)+1,math.floor(img_size[1]/2)+1];
    ring.radius_init = 10;
    ring.halfwidth = 5;
    ring.radius_end = np.amin(ring.center)-1-ring.halfwidth;
    ring.radius = list(range(ring.radius_init, ring.radius_end,2*ring.halfwidth));
    ring.number = len(ring.radius);
    return ring

def generate_p_value_mask (img,ring,options):
    p_value_mask = np.empty(img.shape);
    p_value_mask[:] = np.nan;
    pos_x=[None] * ring.number;
    pos_y=[None] * ring.number;
    angles=[None] * ring.number;
    fft_in_ring=[None] * ring.number;
    outray_ind=[None] * ring.number;
    inray_ind=[None] * ring.number;
    Mu=[None] * ring.number;
    Sigma=[None] * ring.number;
    p_valuesImg=[None] * ring.number;
    for i in range(0,ring.number):
        [pos_x[i],pos_y[i],angles[i]] = get_ring_pos(ring.radius[i],ring.halfwidth,ring.center);
        fft_in_ring[i] = img[pos_x[i],pos_y[i]];
        [outray_ind[i],inray_ind[i]] = get_ray_ind(angles[i],options);
        [Mu[i],Sigma[i],p_valuesImg[i]] = fit_2D_gaussian(fft_in_ring[i],outray_ind[i]);
        # p_value_mask(pos_x[i],pos_y[i]) = p_valuesImg[i];
        
    # statistics=namedtuple('statistics',['Mu','Sigma','p_values']);
    # statistics.Mu = Mu;
    # statistics.Sigma = Sigma;
    # statistics.p_values = p_values;
    
    # positions=namedtuple('positions',['pos_x','pos_y','angles','fft_in_ring']);
    # positions.pos_x = pos_x;
    # positions.pos_y = pos_y;
    # positions.angles = angles;
    # positions.fft_in_ring = fft_in_ring;
    
    # indecies=namedtuple('indecies',['outray_ind','inray_ind','p_values']);
    # indecies.outray_ind = outray_ind;
    # indecies.inray_ind = inray_ind;
    
    return  p_value_mask #, statistics, positions, indecies

def get_ring_pos(ring_radius,ring_halfwidth,center):
    x = list(range(-ring_radius-ring_halfwidth, ring_radius+ring_halfwidth+1,1))
    y = list(range(-ring_radius-ring_halfwidth, ring_radius+ring_halfwidth+1,1))

    [X,Y]=np.meshgrid(x,y);
    XY = np.add(np.power(X,2), np.power(Y,2));
    Norms = np.sqrt(XY);
    trueNorms = (Norms >= [ring_radius-ring_halfwidth]) & (Norms <= [ring_radius+ring_halfwidth]);
    [indX,indY]= np.where(trueNorms)
    pos_x = X[indY,indX];
    pos_y = Y[indY,indX];
    angles = np.arctan2(pos_x,pos_y);
    pos_x = center[0]+X[indY,indX]-1;
    pos_y = center[1]+Y[indY,indX]-1;
    return pos_x, pos_y, angles

def get_ray_ind(angles,options):
    cond1 = (angles >= options.angle_threshold) & (angles <= math.pi - options.angle_threshold);
    cond2 = (angles >= options.angle_threshold - math.pi) & (angles <= -options.angle_threshold);
    cond3 = (angles >= -options.angle_threshold) & (angles <= options.angle_threshold);
    cond4 = angles >= (math.pi - options.angle_threshold);
    cond5 = angles <= (-math.pi + options.angle_threshold);
    
    outray_ind=np.concatenate((np.where(cond1)[0], np.where(cond2)[0]), axis=None);

    inray_ind = np.concatenate((np.where(cond3)[0],np.where(cond4)[0],np.where(cond5)[0]), axis=None);
    
    return outray_ind, inray_ind

def fit_2D_gaussian(fft_in_ring,outray_ind):
    
    siz = fft_in_ring.shape;
    if len(siz)>1:
        if siz[0] < siz[1]: fft_in_ring = np.transpose(fft_in_ring);
        
    siz = outray_ind.shape;
    if len(siz)>1:
        if siz[0] < siz[1]: outray_ind = np.transpose(outray_ind);
	
    Data = [np.real(fft_in_ring),np.imag(fft_in_ring)];
    Data_outray = [np.real(fft_in_ring[outray_ind]),np.imag(fft_in_ring[outray_ind])];
    
    Mu = np.mean(Data_outray, axis=1);#WEIRD, the mean is different to matlab mean
    Sigma = np.cov(Data_outray);
    Data=np.subtract(np.transpose(Data), np.resize(Mu,(len(Data[0]),len(Mu))))
    Data = np.matmul(Data,scipy.linalg.fractional_matrix_power(Sigma,-0.5));
    Data_rayl = np.sum(np.power(Data,2),axis=1);
    p_values = rayleigh.cdf(Data_rayl,q='upper');#we need to get the upper tail probability

    return Mu,Sigma,p_values

def generate_binary_mask(p_value_mask,options):
    
    binary_mask = np.zeros(p_value_mask.shape);
    p_value_mask_filtered = medfilt2d(p_value_mask,kernel_size=options.medfilt_size);
    binary_mask[np.where(p_value_mask_filtered <= options.cutoff_threshold)] = 1;
    
    return binary_mask
