function [] = write_geometryfile (ptcs, filename)
   [interfaces, boundaries] = nrbmultipatch (ptcs);

   % magnetic problem
   boundaries = [];
   boundaries(1).patches = [1 1 2 3 3 4 4 5 6 6];
   boundaries(1).faces   = [2 3 2 2 4 1 4 1 1 3];

   for ibnd = 1:length(boundaries)
      boundaries(ibnd).nsides = length(boundaries(ibnd).patches);
   end

   nrbexport (ptcs, interfaces, boundaries, [filename '_mag.txt']);

   % mechanic problem
   boundaries = [];
   boundaries(1).patches = [5];
   boundaries(1).faces   = [1];
   % what about 3/4
   boundaries(2).patches = [5 5 5];
   boundaries(2).faces   = [2 3 4];

   for ibnd = 1:length(boundaries)
      boundaries(ibnd).nsides = length(boundaries(ibnd).patches);
   end

   nrbexport (ptcs, interfaces, boundaries, [filename '_mech.txt']);
end
