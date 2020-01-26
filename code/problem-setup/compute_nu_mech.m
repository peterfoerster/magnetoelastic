function [nu] = compute_nu_mech (x, y, iptc, nu)
   switch (iptc)
      case {5}
         % plate
         nu = nu * ones(size(x));
      otherwise
         nu = zeros(size(x));
   end%switch
end
