function [] = plotres(img,rec)
if ~isreal(rec)
	warning('showing the real part of the recovered image');
end
fig1 = figure;
ax(1)=subplot(1,2,1); imagesc(img); colormap gray; axis image; title('original image');%caxis([0,255]);
ax(2)=subplot(1,2,2); imagesc(real(rec)); colormap gray; axis image; title('recovered image');%caxis([0,255]);
linkaxes(ax,'xy');
fig2 = figure;
imagesc(abs(img- real(rec))); colormap gray; axis image; title('|img-rec|');
end