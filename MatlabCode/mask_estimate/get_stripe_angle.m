function [angle,binom_val] = get_stripe_angle(p_value_mask, ring, options)

if nargin<3
	options.eps = 1e-3;
	options.drawres = false;
end

size_min = ring.center - floor(ring.radius_end/sqrt(2));
size_max = ring.center + floor(ring.radius_end/sqrt(2));

p_values = p_value_mask(size_min(1):size_max(1),size_min(2):size_max(2));
p_values(isnan(p_values)==1)= 1;

p_values(p_values<=1e-32) = 1e-32;
center = floor(size(p_values)/2)+1;
p_val = p_values;

% imagesc(p_val);

eps = options.eps;
p=length(find(p_val<=eps))/numel(p_val);
[siz1,siz2]=size(p_val);

angle = linspace(-pi/2+0.01,pi/2-0.01,180);
mask=zeros(size(p_val));
if options.drawres
	figure;
	im = imagesc(mask);
	drawnow;
end

for i=1:length(angle)
    if angle(i)>=-pi/4 && angle(i)<0
    X = siz2-center(2);
    Y = center(1)-floor(X*tan(angle(i)));
    [x,y] = bresenham(center(1),center(2),siz1,Y);
    ray = p_val(sub2ind(size(p_val),y,x));
    v = length(find(ray<=eps));
    N = numel(ray);
    binom_val(i) = betainc(p,v+1,N-v);
    if options.drawres
	    mask = zeros(size(p_val));
	    mask(sub2ind(size(mask),y,x))=1;
	    im.CData=mask;
	    drawnow;
	end
    end

    if angle(i)>=0 && angle(i)<pi/4
    X = siz2-center(2);
    Y = center(1)-floor(X*tan(angle(i)));
    [x,y] = bresenham(center(1),center(2),siz1,Y);
    ray = p_val(sub2ind(size(p_val),y,x));
    v = length(find(ray<=eps));
    N = numel(ray);
    binom_val(i) = betainc(p,v+1,N-v);
    if options.drawres
	    mask = zeros(size(p_val));
	    mask(sub2ind(size(mask),y,x))=1;
	    im.CData=mask;
	    drawnow;
	end
    end

    if angle(i)<-pi/4 && angle(i)>=-pi/2
    Y = siz2-center(2);
    X = center(1)-floor(Y*tan(pi/2-angle(i)));
    [x,y] = bresenham(center(1),center(2),X,siz2);
    ray = p_val(sub2ind(size(p_val),y,x));
    v = length(find(ray<=eps));
    N = numel(ray);
    binom_val(i) = betainc(p,v+1,N-v);
    if options.drawres
	    mask = zeros(size(p_val));
	    mask(sub2ind(size(mask),y,x))=1;
	    im.CData=mask;
	    drawnow;
    end
	end

    if angle(i)>pi/4 && angle(i)<=pi/2
    Y = siz2-center(2);
    X = center(1)+floor(Y*tan(pi/2-angle(i)));
    [x,y] = bresenham(center(1),center(2),X,1);
    ray = p_val(sub2ind(size(p_val),y,x));
    v = length(find(ray<=eps));
    N = numel(ray);
    binom_val(i) = betainc(p,v+1,N-v);
    if options.drawres
	    mask = zeros(size(p_val));
	    mask(sub2ind(size(mask),y,x))=1;
	    im.CData=mask;
	    drawnow;
	end
    end
end
if options.drawres
	figure;
	polar(angle,binom_val);
	figure;
	plot(rad2deg(angle),log10(binom_val));
end
end