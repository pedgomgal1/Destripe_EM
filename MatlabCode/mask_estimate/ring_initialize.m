function [ring] = ring_initialize(img_siz)
	ring.center = floor(img_siz/2)+1;
	ring.radius_init = 10;
	ring.halfwidth = 5;
	ring.radius_end = min(ring.center)-1-ring.halfwidth;
	ring.radius = ring.radius_init:2*ring.halfwidth:ring.radius_end;
	ring.number = length(ring.radius);
end