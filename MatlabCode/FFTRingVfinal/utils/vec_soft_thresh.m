function [Xh,Xv] = vec_soft_thresh(Yh,Yv,beta,mode)
    switch lower(mode)
    	case 'l2'
    		modY = sqrt(abs(Yh).^2+abs(Yv).^2);
		    Xh = Yh - min(beta,modY) .* (Yh ./ modY);
		    Xv = Yv - min(beta,modY) .* (Yv ./ modY);
		    
		    Xh(modY == 0) = 0;
		    Xv(modY == 0) = 0;
		case 'l1'
			modYh = abs(Yh);
			modYv = abs(Yv);
			Xh = Yh - min(beta,modYh) .* (Yh ./ modYh);
		    Xv = Yv - min(beta,modYv) .* (Yv ./ modYv);
		    
		    Xh(modYh == 0) = 0;
		    Xv(modYv == 0) = 0;
		case 'linf'
			Xh = Yh;
			Xv = Yv;

			modYh = abs(Yh);
			modYv = abs(Yv);

			ind_Xh = find(modYh>modYv);
			ind_Xv = find(modYh<=modYv);

			Xh(ind_Xh) = Yh(ind_Xh) - min(beta,modYh(ind_Xh)) .* (Yh(ind_Xh)./ modYh(ind_Xh));
			Xv(ind_Xv) = Yv(ind_Xv) - min(beta,modYv(ind_Xv)) .* (Yv(ind_Xv)./ modYv(ind_Xv));
			
			Xh(modYh(ind_Xh) == 0) = 0;
		    Xv(modYv(ind_Xv) == 0) = 0;

    	otherwise
    		
    end
    
end
