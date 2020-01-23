% Orthotropic material.

% The elastic properties of a linear elastic orthotropic material
% are defined by 9 constants: E1, E2, E3, mu12, mu13, mu23, nu12, nu13,
% nu23. (see
% https://en.wikipedia.org/wiki/Orthotropic_material#Stiffness_and_compliance_matrices_in_orthotropic_elasticity)

% We also have to define the three orthogonal directions that determine
% the orthotropic directions.
% I set these directions by calling compute_D and compute_C with the
% material directions for every point (DF, that is 3x3xn_pts array).
% The entries DF(:, i, :) correspond to the Euclidean coordinates (first
% index), of the i-th vector (second index), for all the points (third
% index). I will call these three vectors as a1, a2, a3.

% a1, a2 and a3 are orthonormalized in the following way:
% 
%  b1 = a1 / |a1|
%  b3 = a1 x a2 / |a1 x a2|
%  b2 = a3 x a1

% The directions used for the othotropy (e1,e2,e3) are computed as
% rotation of the vectors b1,b2,b3, aroung the direction b3. I.e.
%
%  e1 =  cos(angle_xy) * b1  +  sin(angle_xy) * b2
%  e2 = -sin(angle_xy) * b1  +  cos(angle_xy) * b2
%  e3 = b3
% 



% Material definition.
material = OrthotropicMaterial;
material.E1   = 25;
material.E2   =  1;
material.E3   =  1;
material.mu12 =  0.5;
material.mu13 =  0.5;
material.mu23 =  0.2;
material.nu13 = 0.25;
material.nu12 = 0.25;
material.nu23 = 0.25;
material.angle_xy = pi/6.0 ;


n_pts = 10;

% Local directions (I usually use here the covariant derivatives of the map).
DF = zeros(3, 3, n_pts);
DF(1, 1, :) = 1.0;
DF(2, 2, :) = 1.0;
DF(3, 3, :) = 1.0;


D = material.compute_D(DF); % 6x6 Voigt matrix at every point.
C = material.compute_C(DF); % 3x3x3x3 elasticity tensor at every point.
