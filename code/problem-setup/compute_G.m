function [G_iptc] = compute_G (x, y, iptc, G)
   switch (iptc)
      case {5}
         % plate
         G_iptc = G * ones(size(x));
      otherwise
         G_iptc = zeros(size(x));
   end%switch
end
