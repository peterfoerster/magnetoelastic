function [mu_ptc] = compute_mu (x, y, iptc, mu)
   switch (iptc)
      case {5}
         % plate
         mu_ptc = mu * ones(size(x));
      otherwise
         % vacuum
         mu_ptc = 4*pi*1e-7 * ones(size(x));
   end%switch
end
