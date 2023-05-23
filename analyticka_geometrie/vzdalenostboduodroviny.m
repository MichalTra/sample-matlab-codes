function [v,d] = vzdalenostboduodroviny(xyz,bod,vect)
% bod a normalovy vect(or) definuji rovinu; je take mozne zadat rovinu
% pomoci rovnice roviny (dat ji jako bod a vect nezadavat/ponechat prazdny)

% v je vektor, jak se od roviny dostat k bodu (xyz), d je vzdalenost

% last modified: 15.3.2019
% category: math

if nargin<2
    error('Not enough input parameters');
end

% vect chci normalizovany
if nargin==2 || isempty(vect)
    vect = bod(1:3)/norm(bod(1:3));
    dd = bod(4)/norm(bod(1:3)); % kdyz normalizuju vect, zmeni se mi i dd
else
    vect = vect/norm(vect);
    dd = -dot(bod,vect); % rovnice roviny je vect(1)*x + vect(2)*y + vect(3)*z + dd = 0
end

v = [xyz ones(size(xyz,1),1)]*[vect(:); dd]/sum(vect.^2);
v = [vect(1)*v vect(2)*v vect(3)*v];

d = sqrt(sum(v.^2,2));

end