function [domain] = discretize_domain (plate, outer_boundary)
   p1 = nrbeval(plate(1), 0);
   p2 = nrbeval(outer_boundary(1), 1);
   domain(1) = nrbline(p1, p2);

   p1 = nrbeval(plate(1), 1);
   p2 = nrbeval(outer_boundary(2), 1);
   domain(2) = nrbline(p1, p2);

   p1 = nrbeval(plate(1), 1);
   p2 = nrbeval(outer_boundary(4), 0);
   domain(3) = nrbline(p1, p2);

   p1 = nrbeval(outer_boundary(8), 1);
   p2 = nrbeval(plate(1), 0);
   domain(4) = nrbline(p1, p2);
end
