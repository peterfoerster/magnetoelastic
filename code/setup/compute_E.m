function [E_iptc] = compute_E (x, y, iptc, E)
   switch (iptc)
      case {5}
         % plate
         E_iptc = E * ones(size(x));
      otherwise
         E_iptc = zeros(size(x));
   end%switch
end
