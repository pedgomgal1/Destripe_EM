function [Mu,Sigma,p_values] = fit_2D_gaussian(fft_in_ring,outray_ind)
	if iscell(fft_in_ring) || iscell(outray_ind)
		warning('the input to fit_2D_gaussian is a cell not a vector');
	end
	siz = size(fft_in_ring);
	if siz(1)<siz(2)
		warning('the first input to fit_2D_gaussian is not column vector.');
		fft_in_ring = transpose(fft_in_ring);
	end
	siz = size(outray_ind);
	if siz(1)<siz(2)
		warning('the second input to fit_2D_gaussian is not column vector.');
		outray_ind = transpose(outray_ind);
	end
	Data = [real(fft_in_ring),imag(fft_in_ring)];
	Data_outray = [real(fft_in_ring(outray_ind)),imag(fft_in_ring(outray_ind))];
	Mu = mean(Data_outray);
	Sigma = cov(Data_outray);
	Data = Data - repmat(Mu,length(Data),1);
	Data = Data * Sigma^(-.5);
	Data_rayl = sum(Data.^2,2);
	p_values = raylcdf(Data_rayl,1,'upper');
end