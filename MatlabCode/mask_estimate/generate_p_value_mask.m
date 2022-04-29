function [p_value_mask, statistics, positions, indecies ] = generate_p_value_mask(img_fft_shift,ring,options)
	p_value_mask = nan(size(img_fft_shift));
	for i = 1:ring.number
		[pos_x{i},pos_y{i},angles{i}] = get_ring_pos(ring.radius(i),ring.halfwidth,ring.center);
		fft_in_ring{i} = img_fft_shift(sub2ind(size(img_fft_shift),pos_x{i},pos_y{i}));
		[outray_ind{i},inray_ind{i}] = get_ray_ind(angles{i},options);
		[Mu{i},Sigma{i},p_values{i}] = fit_2D_gaussian(fft_in_ring{i},outray_ind{i});
		p_value_mask(sub2ind(size(img_fft_shift),pos_x{i},pos_y{i})) = p_values{i};
	end
	statistics.Mu = Mu;
	statistics.Sigma = Sigma;
	statistics.p_values = p_values;

	positions.pos_x = pos_x;
	positions.pos_y = pos_y;
	positions.angles = angles;
	positions.fft_in_ring = fft_in_ring;

	indecies.outray_ind = outray_ind;
	indecies.inray_ind = inray_ind;
end