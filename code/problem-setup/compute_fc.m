function [fc] = compute_fc (x, y, iptc, f)
   switch (iptc)
      case {5}
         % plate
         fc = f * ones(size(x));
      otherwise
         fc = zeros(size(x));
   end%switch
end
