function [nu] = compute_nu (x, y, iptc, nu)
   switch (iptc)
      case {5}
         % anisotropic
         nu = nu * ones(size(x));
      otherwise
         % vacuum
         nu = 1/(pi*4e-7) * ones(size(x));
   end%switch
end
