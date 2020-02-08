function [plate] = create_plate_v2 ()
   p1 = [3 3];
   p2 = [3 4];
   plate(1) = nrbline(p1, p2);

   p1 = [0 4];
   p2 = nrbeval(plate(1), 1);
   plate(2) = nrbline(p1, p2);

   p1 = [0 3];
   p2 = nrbeval(plate(2), 0);
   plate(3) = nrbline(p1, p2);

   p1 = nrbeval(plate(3), 0);
   p2 = nrbeval(plate(1), 0);
   plate(4) = nrbline(p1, p2);
end
