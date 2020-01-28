function [f_iptc] = compute_f (x, y, iptc, f)
   switch (iptc)
      case {5}
         % plate
         f_iptc = f * ones(size(x));
      otherwise
         f_iptc = zeros(size(x));
   end%switch
end
