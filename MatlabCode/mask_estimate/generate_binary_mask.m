function [binary_mask] = generate_binary_mask(p_value_mask,options)
switch options.mask_mode
case 'thresh'
	binary_mask = zeros(size(p_value_mask));
	binary_mask(p_value_mask<=options.cutoff_threshold) = 1;
case 'open'
	binary_mask = zeros(size(p_value_mask));
	binary_mask(p_value_mask<=options.cutoff_threshold) = 1;
	se = strel(options.open_shape,options.open_size);
	binary_mask = imopen(binary_mask,se);
case 'medfilt'
	binary_mask = zeros(size(p_value_mask));
	p_value_mask_filtered = medfilt2(p_value_mask,options.medfilt_size);
	binary_mask(p_value_mask_filtered<=options.cutoff_threshold) = 1;
end
end