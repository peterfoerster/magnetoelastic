function [kappa_ptc] = compute_kappa (x, y, iptc, omega, sigma, epsilon)
   i_cmp = complex(0, 1);
   switch (iptc)
      case {5}
         % plate
         kappa_ptc = (i_cmp*omega*sigma - omega^2 * epsilon) * ones(size(x));
      otherwise
         kappa_ptc = zeros(size(x));
   end%switch
end
