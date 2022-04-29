function [F,outlier_ind] = replace_outiler(img_fft_shift, binary_mask, statistics, positions, indecies, ring, options)

	F = img_fft_shift;
	outlier_ind = find(binary_mask==1);
	[outlier_x,outlier_y] = ind2sub(size(img_fft_shift),outlier_ind);
	X = outlier_x - ring.center(1);
	Y = outlier_y - ring.center(2);
	Norms = sqrt(X.^2+Y.^2);
	ring_up = ring.radius+ring.halfwidth;
	ring_down = ring.radius-ring.halfwidth;

	for i = 1:length(outlier_x)/2
		ring_index = find(ring_up>Norms(i)&ring_down<=Norms(i));
        if isempty(ring_index)
           ring_index = find(ring_up>=Norms(i)&ring_down<Norms(i)); 
        end
		if options.random_on
			temp = mvnrnd(statistics.Mu{ring_index},statistics.Sigma{ring_index});
			sample = temp(1)+1j*temp(2);
			clear temp;
		else
			sample = sqrt(statistics.Mu{ring_index}(1)^2+statistics.Sigma{ring_index}(1,1)+...
						  statistics.Mu{ring_index}(2)^2+statistics.Sigma{ring_index}(2,2));
		end
		F(outlier_x(i),outlier_y(i)) = sample;
		F(outlier_x(end-i+1),outlier_y(end-i+1)) = conj(sample);
	end

end