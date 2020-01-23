function [outer_boundary] = create_outerboundary (plate)
   p1 = [5 0];
   p2 = nrbeval(plate(1), 0);
   p2 = [5 p2(2)];
   outer_boundary(1) = nrbline(p1, p2);

   p1 = nrbeval(outer_boundary(1), 1);
   p2 = nrbeval(plate(1), 1);
   p2 = [p1(1) p2(2)];
   outer_boundary(2) = nrbline(p1, p2);

   p1 = nrbeval(outer_boundary(2), 1);
   p2 = [5 6];
   outer_boundary(3) = nrbline(p1, p2);

   p1 = nrbeval(plate(1), 1);
   p2 = nrbeval(outer_boundary(3), 1);
   p1 = [p1(1) p2(2)];
   outer_boundary(4) = nrbline(p1, p2);

   p1 = [0 6];
   p2 = nrbeval(outer_boundary(4), 0);
   outer_boundary(5) = nrbline(p1, p2);

   p1 = nrbeval(plate(3), 1);
   p2 = nrbeval(outer_boundary(5), 0);
   outer_boundary(6) = nrbline(p1, p2);

   p1 = [0 0];
   p2 = nrbeval(plate(3), 0);
   outer_boundary(7) = nrbline(p1, p2);

   p1 = nrbeval(outer_boundary(7), 0);
   p2 = nrbeval(plate(4), 1);
   p2 = [p2(1) 0];
   outer_boundary(8) = nrbline(p1, p2);

   p1 = nrbeval(outer_boundary(8), 1);
   p2 = nrbeval(outer_boundary(1), 0);
   outer_boundary(9) = nrbline(p1, p2);
end
