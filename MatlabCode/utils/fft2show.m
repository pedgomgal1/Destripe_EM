function fft2show(img)
fig = figure;
imagesc(log10(abs(fftshift(fft2(img)))));
colormap('gray');
colorbar;
% figure;
% imagesc(angle(fftshift(fft2(img))));
% colormap('gray');
% colorbar;
end

