function [tau] = compute_tau (x, y, ib, tau)
   switch (ib)
      case {2}
         % plate
         tau = [tau * ones(size(x)); zeros(size(x))];
         % tau = [zeros(size(x)); tau * ones(size(x))];
      otherwise
         tau = [zeros(size(x)); zeros(size(x))];
   end%switch
end
