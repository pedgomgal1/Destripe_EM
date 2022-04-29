function [Z ] = FFT2_TV3_L_Proj(img,mask_fft,options)

% This is the main solver to be used for the FFTRingV2.
% inputs:
% 		img: 3D image Volume
%		mask_fft: 3D volume of fft-mask for each slice
%		options: options for the optimizer. can leave empty
% 		for default values.
% Output: 
% 		Z: recovered image (Denoised)
% Note: Addaptive step size slowes the code down a bit because
% of extra calculations but should have better convergence properties.

% Default values for options
if nargin<3
	fprintf('Default initialization of parameters are set!\n');
	rho_Z = 1;
	rho_X = 1;
	rho_Z_tild = 1;
	scale = 1e7;
	lambda_X = 1e-2;
	lambda_Ztild = 1e-1;
	threshold_mode = 'l2';
	maxiter = 100;
	tol = 1e-4;
	verbose = true;
	adaptive = true;
else
	if isfield(options,'rho_Z')
		rho_Z = options.rho_Z;
	else
		rho_Z = 1;
	end
	if isfield(options,'rho_X')
		rho_X = options.rho_X;
	else
		rho_X = 1;
	end
	if isfield(options,'rho_Z_tild')
		rho_Z_tild = options.rho_Z_tild;
	else
		rho_Z_tild = 1;
	end
	if isfield(options,'scale')
		scale = options.scale;
	else
		scale = 1e7;
	end
	if isfield(options,'lambda_X')
		lambda_X = options.lambda_X;
	else
		lambda_X = 1e-2;
	end
	if isfield(options,'lambda_Ztild')
		lambda_Ztild = options.lambda_Ztild;
	else
		lambda_Ztild = 1e-1;
	end
	if isfield(options,'threshold_mode')
		threshold_mode = options.threshold_mode;
	else
		threshold_mode = 'l2';
	end
	if isfield(options,'maxiter')
		maxiter = options.maxiter;
	else
		maxiter = 100;
	end
	if isfield(options,'tol')
		tol = options.tol;
	else
		tol = 1e-4;
	end
	if isfield(options,'verbose')
		verbose = options.verbose;
	else
		verbose = true;
	end
	if isfield(options,'adaptive')
		adaptive = options.adaptive;
	else
		adaptive = false;
	end
end


W = double(~mask_fft) .* (ones(size(img)));
clear('mask_fft'); % clearing to get memory efficiency

% setup
fftw('planner','measure');
Y = img;
Y_hat = fft2(Y);
upper_thresh = 1;
lower_thresh = 0;
clear('img');

% define the difference operator kernel
Dh = @(X) X - [X(:,end,:),X(:,1:end-1,:)];
DhT = @(X) X - [X(:,2:end,:),X(:,1,:)];

Dv = @(X) X - [X(end,:,:);X(1:end-1,:,:)];
DvT = @(X) X - [X(2:end,:,:);X(1,:,:)];

Dd = @(X) X - cat(3,X(:,:,end),X(:,:,1:end-1));
DdT = @(X) X - cat(3,X(:,:,2:end),X(:,:,1));

% define the difference operator kernel
dh = zeros(size(Y));
dh(1,1,1) = 1;
dh(1,2,1) = -1;

dv = zeros(size(Y));
dv(1,1,1) = 1;
dv(2,1,1) = -1;

dd = zeros(size(Y));
dd(1,1,1) = 1;
dd(1,1,2) = -1;

lv = zeros(size(Y));
lv(1,1,1)=1;
lv(2,1,1)=-2;
lv(3,1,1)=1;

% compute FFTs for filtering
FDH = fftn(dh);
FDV = fftn(dv);
FDD = fftn(dd);
F2D = abs(FDH).^2 + abs(FDV).^2 + abs(FDD).^2;
FLV = fftn(lv);
F2LV = abs(FLV).^2;
F2LV_FY = F2LV.*fftn(Y);
clear('FDH', 'FDV', 'FDD', 'FLV','dh','dv','dd','lv'); 


% Initialization
Z = zeros(size(Y)); % we can initialize with zeros!
Z_tild = Z;
Z_hat = fft2(Z)/scale;
Xh = Dh(Z_tild);
Xv = Dv(Z_tild);
Xd = Dd(Z_tild);

% Dual variable initializations
Uz = zeros(size(Z));
Utild = zeros(size(Z_tild));
UXh = zeros(size(Xh));
UXv = zeros(size(Xv));
UXd = zeros(size(Xd));

% adaptive step size parameters
alpha_adapt = 0.7;
gamma_adapt = 2;

clear('Y');

% initialize the gaps -- 
% gaps are used for adaptive step and convergence (Duality Gap should be small at optimal)
gapZ_Zhat(1) = 0;
gapZ_Ztild(1) = 0;
gapX_Ztild(1) = 0;

if verbose
	obj1(1) = sum(sum(sum((W.*(abs((scale*Z_hat)-Y_hat)).^2))));
	switch lower(threshold_mode)
		case 'l2'
			obj2(1) = sum(sum(sum(sqrt(Xh.^2+Xv.^2+Xd.^2))));	
		case 'l1'
			obj2(1) = sum(sum(sum(abs(Xh)+abs(Xv)+abs(Xd))));
		case 'linf'
			obj2(1) = sum(sum(sum(max(abs(Xh),abs(Xv),abs(Xd)))));
		otherwise
			error('threshold mode not defined');
	end
	fprintf('i:%4d|o1:%2.2e o2:%2.2e|G1:%2.2e G2:%2.2e G3:%2.2e|\n',0,obj1(1),obj2(1),gapZ_Zhat(1),gapZ_Ztild(1),gapX_Ztild(1));
end

tic
for i = 1:maxiter
	% Update Z_hat ---> (scale diag(Omega)+rho_Z/2I)^-1(diag(Omega).*Y_hat+rho_Z/2 F(Z-U)) diag solve in fft2 domain
	Z_hat = (W.*Y_hat + (rho_Z/2) * fft2(Z+Uz)) ./ (scale*W + (rho_Z/2)*scale);

	Z_ifft  = real(ifft2(Z_hat*scale));
	% Update Z_tild ---> diag-solve in fftn domian
	Z_tild = ifftn((lambda_Ztild*F2LV_FY + fftn( (rho_X/2)*DhT(Xh+UXh) + (rho_X/2)*DvT(Xv+UXv) + (rho_X/2)*DdT(Xd+UXd) + (rho_Z_tild/2)*(Z+Utild)))...
			 ./ (lambda_Ztild*F2LV +(rho_X/2)*F2D + (rho_Z_tild/2)));

	% Update Xh, Xv ----> vector soft-thresholding
	NXh     = Dh(Z_tild)-UXh;
	NXv     = Dv(Z_tild)-UXv;
	NXd     = Dd(Z_tild)-UXd;
    [Xh,Xv,Xd] = vec_soft_thresh3D(NXh,NXv,NXd, lambda_X/rho_X,threshold_mode);

    % Projection update
    Z = (rho_Z*(Z_ifft - Uz)+rho_Z_tild*(Z_tild-Utild))/(rho_Z+rho_Z_tild);
    Z(Z<lower_thresh) = lower_thresh;
    Z(Z>upper_thresh) = upper_thresh;

    % dual update
    DhZ_tild_saved = Dh(Z_tild);
    DvZ_tild_saved = Dv(Z_tild);
    DdZ_tild_saved = Dd(Z_tild);
    UXh = UXh + (-DhZ_tild_saved+Xh);
	UXv = UXv + (-DvZ_tild_saved+Xv);
	UXd = UXd + (-DdZ_tild_saved+Xd);
	Uz = Uz + (-Z_ifft+Z);
	Utild = Utild +(-Z_tild+Z);

	% Update Gap values
	gapZ_Zhat(i+1) = norm(Z_ifft(:)-Z_hat(:),'fro');
	gapZ_Ztild(i+1) = norm(Z(:)-Z_tild(:),'fro');
	gapX_Ztild(i+1) = norm((-DhZ_tild_saved(:)+Xh(:)).^2+(-DvZ_tild_saved(:)+Xv(:)).^2+(-DdZ_tild_saved(:)+Xd(:)).^2,'fro');

	% Update rho values w.r.t. the gap if adaptive is required
	if adaptive
		if gapZ_Zhat(i+1)>= alpha_adapt*gapZ_Zhat(i);
    		rho_Z = gamma_adapt*rho_Z;
    	end
    	if gapZ_Ztild(i+1)>= alpha_adapt*gapZ_Ztild(i);
    		rho_Z_tild = gamma_adapt*rho_Z_tild;
    	end
    	if gapX_Ztild(i+1)>= alpha_adapt*gapX_Ztild(i);
    		rho_X = gamma_adapt*rho_X;
    	end
	end

	if verbose
		obj1(i+1) = sum(sum(sum((W.*(abs((scale*Z_hat)-Y_hat)).^2))));
		switch lower(threshold_mode)
			case 'l2'
				obj2(i+1) = sum(sum(sum(sqrt(Xh.^2+Xv.^2+Xd.^2))));	
			case 'l1'
				obj2(i+1) = sum(sum(sum(abs(Xh)+abs(Xv)+abs(Xd))));
			case 'linf'
				obj2(i+1) = sum(sum(sum(max(abs(Xh),abs(Xv),abs(Xd)))));
			otherwise
				error('threshold mode not defined');
		end
		fprintf('i:%3d|o1:%2.2e o2:%2.2e|G1:%2.2e G2:%2.2e G3:%2.2e|\n',i,obj1(i+1),obj2(i+1),gapZ_Zhat(i+1),gapZ_Ztild(i+1),gapX_Ztild(i+1));
	end 

	if (gapZ_Ztild(i+1)<tol) && (gapX_Ztild(i+1)<tol) && (abs(gapZ_Zhat(i+1)-gapZ_Zhat(i))<tol)
		break;
	end

end

time_all = toc;

% summary of the results 
fprintf('--------------------------------------------------------------\n');
fprintf('Summary of Results\n');
fprintf('--------------------------------------------------------------\n');
R1 = abs(Z-Z_tild);
R2h = abs(Xh - Dh(Z_tild));
R2v = abs(Xv - Dv(Z_tild));
R2d = abs(Xd - Dd(Z_tild));
R3 = abs(Z - real(ifft2(Z_hat*scale)));
fprintf('max{|Z - Z_tild|} = %4.4e\n',max(R1(:)));
fprintf('max{|Xh - Dh(F^-1 Z_hat)|} = %4.4e\n',max(R2h(:)));
fprintf('max{|Xv - Dv(F^-1 Z_hat)|} = %4.4e\n',max(R2v(:)));
fprintf('max{|Xv - Dd(F^-1 Z_hat)|} = %4.4e\n',max(R2d(:)));
fprintf('max{|Z - F^-1 Z_hat|} = %4.4e\n',max(R3(:)));
fprintf('total time = %4.4f\n',time_all);
fprintf('--------------------------------------------------------------\n');

 
