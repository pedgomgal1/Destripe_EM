function [Z,Z_hat,obj1,obj2,obj3,obj4,time_all]=FFT_TV_L_proj_Solver(img,mask_fft,lambda_X,lambda_Ztild,threshold_mode,maxiter,verbose)

img = double(img)/255;
W = double(~mask_fft) .* (ones(size(img)));
lower_bound = min(img(:));
upper_bound = max(img(:));
% setup
fftw('planner','measure');
Y = img;
Y_hat = fft2(img);

% define the difference operator kernel
Dh = @(X) X - [X(:,end),X(:,1:end-1)];
DhT = @(X) X - [X(:,2:end),X(:,1)];

Dv = @(X) X - [X(end,:);X(1:end-1,:)];
DvT = @(X) X - [X(2:end,:);X(1,:)];


% define the difference operator kernel
dh = zeros(size(img));
dh(1,1) = 1;
dh(1,2) = -1;

dv = zeros(size(img));
dv(1,1) = 1;
dv(2,1) = -1;

lv = zeros(size(img));
lv(1,1)=1;
lv(2,1)=-2;
lv(3,1)=1;

% compute FFTs for filtering
FDH = fft2(dh);
FDHC = conj(FDH);
F2DH = abs(FDH).^2;
FDV = fft2(dv);
FDVC = conj(FDV);
F2DV = abs(FDV).^2;
F2D = F2DH + F2DV;

FLV = fft2(lv);
FLVH = conj(FLV);
F2LV = abs(FLV).^2;

% Parameter setup
rho_Z = .1;
rho_X = .1;
rho_Z_tild = .1;

scale = 1e7;

% Initialization
Z = zeros(size(img));
Z_hat = fft2(Z)/scale;
Xh = Dh(Z);
Xv = Dv(Z);
Z_tild = Z;

Uz = zeros(size(Z));
UXh = zeros(size(Xh));
UXv = zeros(size(Xv));
Utild = zeros(size(Z_tild));

FLTLY = FLVH.*FLV.*Y_hat;

if verbose
	obj1 = sum(sum(W.*(abs((scale*Z_hat)-Y_hat)).^2));
	switch lower(threshold_mode)
		case 'l2'
			obj2 = sum(sum(sqrt(Xh.^2+Xv.^2)));		
		case 'l1'
			obj2 = sum(sum(abs(Xh)+abs(Xv)));
		case 'linf'
			obj2 = sum(sum(max(abs(Xh),abs(Xv))));
		otherwise
			error('threshold mode not defined');
	end
	
	obj3 = sum(sum(min(abs(Z(find(Z<0|Z>1))),abs(Z(find(Z<0|Z>1))-1))));
	obj4 = sum(sum(Dv(Dv(Z_tild-Y)).^2));
	fprintf('i: %4d | o1: %4.4e o2: %4.4e o3: %4.4e o4:%4.4e|\n',0, obj1, obj2,obj3,obj4);
end
tic
for i = 1:maxiter
	% Update Z_hat ----> L2-L2 proplem
	Z_hat = (W.*Y_hat + fft2((rho_Z/2)*(Z+Uz)+ (rho_X/2)*DhT(Xh+UXh) + (rho_X/2)*DvT(Xv+UXv)))./...
			( scale*W + scale*(rho_X/2)*F2D + scale * (rho_Z/2));

	% Update Z_tild
	Z_tild = ifft2((lambda_Ztild*FLTLY+(rho_Z_tild/2)*fft2(Z+Utild))./(rho_Z_tild/2+lambda_Ztild*F2LV));

	Z_ifft  = ifft2(Z_hat*scale);
	% Update Xh, Xv ----> vector soft-thresholding
	NXh     = Dh(Z_ifft)-UXh;
	NXv     = Dv(Z_ifft)-UXv;
    [Xh,Xv] = vec_soft_thresh(NXh,NXv, lambda_X/rho_X,threshold_mode);

    % Projection update
    Z = (rho_Z*(Z_ifft - Uz)+rho_Z_tild*(Z_tild-Utild))/(rho_Z+rho_Z_tild);
    Z(Z<lower_bound) = lower_bound;
    Z(Z>upper_bound) = upper_bound;

    UXh = UXh + (-Dh(Z_ifft)+Xh);
	UXv = UXv + (-Dv(Z_ifft)+Xv);
	Uz = Uz + (-Z_ifft+Z);
	Utild = Utild +(-Z_tild+Z);

	if verbose
	obj1(i) = sum(sum(W.*(abs((scale*Z_hat)-Y_hat)).^2));
	switch lower(threshold_mode)
		case 'l2'
			obj2(i) = sum(sum(sqrt(Xh.^2+Xv.^2)));		
		case 'l1'
			obj2(i) = sum(sum(abs(Xh)+abs(Xv)));
		case 'linf'
			obj2(i) = sum(sum(max(abs(Xh),abs(Xv))));
		otherwise
			error('threshold mode not defined');
	end
	obj3(i) = sum(sum(min(abs(Z_ifft(find(Z_ifft<0|Z_ifft>1))),abs(Z_ifft(find(Z_ifft<0|Z_ifft>1))-1))));
	obj4(i) = sum(sum(Dv(Dv(Z_tild-Y)).^2));
	fprintf('i: %4d | o1: %4.4e o2: %4.4e o3: %4.4e o4: %4.4e |\n',i, obj1(i), obj2(i),obj3(i),obj4(i));
	end
end
time_all = toc;

% summary of the results 
if verbose
fprintf('--------------------------------------------------------------\n');
fprintf('Summary of Results\n');
fprintf('--------------------------------------------------------------\n');
R1 = abs(Z-Z_tild);
R2h = abs(Xh - Dh(ifft2(Z_hat*scale)));
R2v = abs(Xv - Dv(ifft2(Z_hat*scale)));
R3 = abs(Z - ifft2(Z_hat*scale));
fprintf('max{|Z - Z_tild|} = %4.4e\n',max(R1(:)));
fprintf('max{|Xh - Dh(F^-1 Z_hat)|} = %4.4e\n',max(R2h(:)));
fprintf('max{|Xv - Dv(F^-1 Z_hat)|} = %4.4e\n',max(R2v(:)));
fprintf('max{|Z - F^-1 Z_hat|} = %4.4e\n',max(R3(:)));
fprintf('total time = %4.4f\n',time_all);
fprintf('--------------------------------------------------------------\n'); 
end
end