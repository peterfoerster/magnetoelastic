function [E] = compute_E (x, y, iptc, E)
   switch (iptc)
      case {5}
         % plate
         E = E * ones(size(x));
      otherwise
         E = zeros(size(x));
   end%switch
end
