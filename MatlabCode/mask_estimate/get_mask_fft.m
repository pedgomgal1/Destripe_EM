function [mask_fft,p_value_mask, statistics, positions, indecies, ring] = get_mask_fft(img,options)

%% fft calculation
[m,n] = size(img);
im_fft = fft2(img);
img_fft_shift = fftshift(im_fft);

%% Ring parameters
[ring] = ring_initialize(size(img_fft_shift));

%% main body
[p_value_mask, statistics, positions, indecies] = generate_p_value_mask(img_fft_shift,ring,options);
binary_mask_medfilt = generate_binary_mask(p_value_mask,options);
mask_fft = ifftshift(binary_mask_medfilt);

end