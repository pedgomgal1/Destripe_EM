function [pos_x,pos_y,angles] = get_ring_pos(ring_radius,ring_halfwidth,center)
	x =  [-ring_radius-ring_halfwidth:ring_radius+ring_halfwidth];
	y =  [-ring_radius-ring_halfwidth:ring_radius+ring_halfwidth];
	[X,Y] = meshgrid(x,y);
	Norms = sqrt(X.^2+Y.^2);
	ind = find(Norms>=ring_radius-ring_halfwidth & Norms<=ring_radius+ring_halfwidth);
	pos_x = X(ind);
	pos_y = Y(ind);
	pos = [pos_x,pos_y];
	angles = atan2(pos_x,pos_y);
	pos_x = center(1)+X(ind);
	pos_y = center(2)+Y(ind);
end