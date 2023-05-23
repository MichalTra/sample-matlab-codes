function [v,d,abc,a,b] = vzdalenostboduodprimky(xyz,bod,vect)
% bod a vect(or) definuji primku

% v je vektor, jak se z primky dostat k bodu (xyz), d je vzdalenost

% last modified: 10.12.2018
% cateogry: math

% xyz = xyz - repmat(bod,size(xyz,1),1); % pomale
xyz = [xyz(:,1)-bod(1) xyz(:,2)-bod(2) xyz(:,3)-bod(3)];

vect = vect/norm(vect);

while true
    a = rand(1,3);
    % p = dot(a,vect); % prumet do vect % pomale
    p = a*vect(:); % prumet do vect
    a = a-p*vect;
    if norm(a)<0.01 % nez riskovat chyby ze zaokrouhlovani, radeji zkusim jiny vektor
        continue
    else
        a = a/norm(a);
        break;
    end
end

b = cross(a,vect); % je normalizovany, protoze ostatni uz normalizovane jsou
% a a b jsou nejake (vice ci mene arbitrary) vektory, ktere spolecne s
% vektorem vect tvori bazi 3D prostoru

abc = xyz/[a;b;vect];

v = abc(:,1:2) * [a;b];
d = sqrt(sum(abc(:,1:2).^2,2));

end