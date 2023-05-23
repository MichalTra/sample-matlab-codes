function v = samplespheresurf(n,r)
% homogeneous sampling of a sphere surface

% INPUT:
% n ... number of vectors
% r ... radius

% last modified: 12.4.2018
% category: math

if nargin==1
    r = 1; % polomer
end

az = rand(n,1)*2*pi; % azimut

el = rand(n,1); % elevace
el = sign(rand(n,1)-0.5).*asin(el); % random sign multiplied by an appropriate distribution (that's why - is not necessary; 0.5 place points on the equator)


[x,y,z] = sph2cart(az,el,r*ones(n,1));
v = [x y z];