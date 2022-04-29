function Visual_polar(img_fft_shift, fft_modified, statistics, positions, indecies, ring)

ring_position = 1;

for i = 1:ring.number
	new_fft_in_ring{i} = fft_modified(sub2ind(size(fft_modified),positions.pos_x{i},positions.pos_y{i}));
end

fig1 = figure();
ax1 = polaraxes();
subplot(1,2,1,ax1);
p1 = polarscatter(ax1,0,0,'b.');
ax2 = polaraxes();
subplot(1,2,2,ax2);
p2 = polarscatter(ax2,0,0,'b.');

p1.RData = abs(positions.fft_in_ring{ring_position});
p1.ThetaData = positions.angles{ring_position};

p2.RData = abs(new_fft_in_ring{ring_position});
p2.ThetaData = positions.angles{ring_position};

drawnow;

txt = uicontrol(fig1,...
		'Style','text',...
    	'String','ring 1',...
    	'Units', 'Normalized',...
        'Position', [0.2 0.05 0.6 0.05]);
sld = uicontrol(fig1,...
		'Style', 'slider',...
        'Min',1,'Max',ring.number,'Value',1,...
        'Units', 'Normalized',...
        'Position', [0.2 0.0 0.6 0.05],...
        'Callback', {@fun1, img_fft_shift, fft_modified, statistics, positions, indecies, ring, new_fft_in_ring, txt, p1, p2 },...
        'SliderStep', [1/(ring.number-1) 1]);

end

function fun1(hObject, eventdata, img_fft_shift, fft_modified, statistics, positions, indecies, ring, new_fft_in_ring, txt, p1, p2 )
ring_position = floor(get(hObject,'Value'));
p1.RData = abs(positions.fft_in_ring{ring_position});
p1.ThetaData = positions.angles{ring_position};

p2.RData = abs(new_fft_in_ring{ring_position});
p2.ThetaData = positions.angles{ring_position};
txt.String = ['ring ',num2str(ring_position)];
end
