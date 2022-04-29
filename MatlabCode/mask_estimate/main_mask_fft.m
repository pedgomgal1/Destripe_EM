function [mask_fft,p_value_mask, statistics, positions, indecies, ring] = main_mask_fft(img,options);

if nargin<2
	options.allout = 0; 	% consider all outliers
	options.angle_threshold = 0.1;
	options.random_on = 1;
	options.cutoff_threshold = 1e-3;
	options.dialation_mode = 'square';
	options.dialation_size = 5;
	options.mask_mode = 'medfilt';
	options.medfilt_size = [5,5];
	options.box_mask = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%% main body %%%%%%%%%%%%%%%%%%%%%%%%%
if options.box_mask==true
	%% Generate a box as mask! -no auto detection
	center = floor(size(img)/2)+1;
	width = 4;
	box_size = 400;
	middle_skip = 10;
	
	mask_fft = zeros(size(img));
	mask_fft(center(1)-width:center(1)+width,[center(2)-(box_size+middle_skip):center(2)-...
	middle_skip,center(2)+middle_skip:center(2)+box_size+middle_skip])=1;
	mask_fft = ifftshift(mask_fft);
	mask_fft(1,1) = 0;
	return;
end
if numel(size(size(img)))==2;
	[mask_fft,p_value_mask, statistics, positions, indecies, ring] = get_mask_fft(img,options);
	[angle,binom_val] = get_stripe_angle(p_value_mask, ring);
	[min_ind,min_val,binom_val_box] = get_stripe_best_box(p_value_mask,ring);
	d = min_ind(1);
	mask_fft_modif = ifftshift(mask_fft);
	mask_fft_modif(1:ring.center(1)-d-10,:)=0;
	mask_fft_modif(ring.center(1)+d+10:end,:)=0;
	switch lower(options.dialation_mode)
		case 'line'
			SE = strel('line',options.dialation_size,angle);
			mask_fft_modif =  imdilate(mask_fft_modif,SE);
		case 'square'
			SE = strel('square',options.dialation_size);
			mask_fft_modif =  imdilate(mask_fft_modif,SE);
		otherwise
			warning('dialation mode was incorrect! default square was done');
			SE = strel('square',5);
			mask_fft_modif =  imdilate(mask_fft_modif,SE);
	end
	
	%%%
	center = floor(size(img)/2)+1;
	center_box_size = 4;
	mask_fft_modif(center(1)-center_box_size:center(1)+center_box_size,center(2)-center_box_size:center(2)+center_box_size)=0;
	%%%
	mask_fft_modif = fftshift(mask_fft_modif);
	mask_fft = mask_fft_modif;
else
	error('Image input for this function must be 2D image');
end

end