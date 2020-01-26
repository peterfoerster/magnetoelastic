function [G] = compute_G (x, y, iptc, G)
   switch (iptc)
      case {5}
         % plate
         G = G * ones(size(x));
      otherwise
         G = zeros(size(x));
   end%switch
end
