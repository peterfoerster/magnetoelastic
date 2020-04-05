function [nu_iptc] = compute_nu (x, y, iptc, nu)
   switch (iptc)
      case {5}
         % plate
         nu_iptc = nu * ones(size(x));
      otherwise
         nu_iptc = zeros(size(x));
   end%switch
end
