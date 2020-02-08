function [] = write_geometryfile_v2 (ptcs, filename)
   [interfaces, boundaries] = nrbmultipatch (ptcs);

   % magnetic problem
   boundaries = [];
   % homogeneous Dirichlet
   boundaries(1).patches = [1 1 2 3 3 4 6 7 7 8 9 9];
   boundaries(1).faces   = [2 3 2 2 4 4 3 1 3 1 1 4];

   for ibnd = 1:length(boundaries)
      boundaries(ibnd).nsides = length(boundaries(ibnd).patches);
   end

   nrbexport (ptcs, interfaces, boundaries, [filename '_mag.txt']);

   % mechanic problem
   boundaries = [];
   % homogeneous Dirichlet
   boundaries(1).patches = [5];
   boundaries(1).faces   = [1];
   % inhomogeneous Neumann
   boundaries(2).patches = [5];
   boundaries(2).faces   = [2];
   % homogeneous Neumann
   boundaries(3).patches = [5 5];
   boundaries(3).faces   = [3 4];

   for ibnd = 1:length(boundaries)
      boundaries(ibnd).nsides = length(boundaries(ibnd).patches);
   end

   nrbexport (ptcs, interfaces, boundaries, [filename '_mec.txt']);
end
