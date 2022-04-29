function [min_ind,min_val,binom_val_box] = get_stripe_best_box(p_value_mask,ring)

size_min = ring.center - floor(ring.radius_end/sqrt(2));
size_max = ring.center + floor(ring.radius_end/sqrt(2));

p_values = p_value_mask(size_min(1):size_max(1),size_min(2):size_max(2));
p_values(isnan(p_values)==1)= 1;

p_values(p_values<=1e-32) = 1e-32;
center = floor(size(p_values)/2)+1;
p_val = p_values;

w_siz = floor(min(size(p_values))/10);
skip = 10;
eps = 1e-3;
p=length(find(p_val<=eps))/numel(p_val);
[siz1,siz2]=size(p_val);

min_val = 1e5;
min_ind = 0;

for d=1:floor(size(p_values)/20)
   ind=0;
   for pos = 1:skip:siz2-w_siz
       ind=ind+1;
       x = p_val(center(1)-d:center(1)+d,pos:pos+w_siz);
       v = length(find(x<=eps));
       N = numel(x);
       binom_val_box(d,ind) = betainc(p,v+1,N-v);
       if binom_val_box(d,ind)<=min_val
          min_ind = [d,pos];
          min_val = binom_val_box(d,ind);
        end
   end
end


end