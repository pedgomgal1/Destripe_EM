function [im_rec] = FFT_Zero_proj(img,mask_fft,options)

% This is the solver for the simple projection method FFT-RingV1.
% The idea is to set the corrupted region in the frequency
% domain to zero. 
% inputs:
%       img: 3D image Volume
%       mask_fft: 3D volume of fft-mask for each slice
%       The region to be removed should be 1 and outside 0.
%       options: struct containing options for the algorithm
%       can skip options and use default values 
% outputs:
%       im_rec: recovered 3D volume image

if nargin<3
    verbose = true;
    maxiter = 200;
    tol = 1e-10;
else 
    if isfield(options,'verbose')
        verbose = options.verbose;
    else
        verbose = true;
    end
    if isfield(options,'maxiter')
        maxiter = options.maxiter;
    else
        maxiter = 20;
    end
    if isfield(options,'tol')
        tol = options.tol;
    else
        tol = 1e-10;
    end
end

% initialization of the variables
Y = img;
N = size(Y,1);
% H,L are limits for the projection if image is normalized to [0,1]
% can also choose L=0, H=1. 
L=min(Y(:)); 
H=max(Y(:));
Z_hat = fft2(Y); % 2D fourier transform of each slice
Y_hat = Y;

P = zeros(size(Y));
Q = zeros(size(Y));


Z_hat_old = zeros(size(Y));
Y_hat_old = zeros(size(Y));

C=0; % The value set in the frequency box 

% Dykstra's iterations
for i = 1:maxiter
    
    temp_Z = Z_hat + P;
    temp_Z(mask_fft==1)=C;
    
    P = Z_hat + P - temp_Z;
    
    Y_hat = ifft2(temp_Z);
    
    temp_Y = Y_hat +Q;
    temp_Y(temp_Y<=L)=L;
    temp_Y(temp_Y>=H)=H;
    
    Z_hat = fft2(temp_Y);
    Q = Y_hat + Q - temp_Y;
    
    Z_change = norm(Z_hat-Z_hat_old,'fro');
    Y_change = norm(Y_hat-Y_hat_old,'fro');
    if verbose
        fprintf('iter: %3d | Z_hat change: %1.5e | Y change: %1.5e\n',i,Z_change,Y_change);
    end
    Z_hat_old = Z_hat;
    Y_hat_old = Y_hat;
    if Z_change<tol && Y_change<tol
        break;
    end
end

im_rec = Y_hat;

end