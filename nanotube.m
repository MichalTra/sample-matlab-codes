function [g,l] = nanotube(type,d,L)
% will make a nanotube with diameter close to d and length close to L
% d, L in Angstroms; aligned in x axis
% type nanotube/zigzag/armchair

% last modified: 1.10.2015
% category: chemistry

xy = [1.2124 2.8000; 0 2.1000; 0 0.7000; 1.2124 0];
l = [2.4249    4.2000];

o = pi*d; % circumfence

switch(type)
    case 'nanotube'
        f = fix(o./l);
        f = (floor(f)+1).*l;
        if f(2)>f(1)
            l([1 2]) = l([2 1]);
            xy(:,[1 2]) = xy(:,[2 1]);
        end
    case 'zigzag'
        l([1 2]) = l([2 1]);
        xy(:,[1 2]) = xy(:,[2 1]);
    case 'armchair'
        % do nothing
    otherwise
        error('The only types allowed are nanotube/zigzag/armchair');
end

m = floor(L/l(1));
n = floor(o/l(2));

[X,Y] = meshgrid(0:1:m,0:1:n);
X = X(:); Y = Y(:);
X = repmat(X,1,size(xy,1)); X = X'; X = X(:);
Y = repmat(Y,1,size(xy,1)); Y = Y'; Y = Y(:);

g = repmat(xy,length(X)/size(xy,1),1);
g(:,1) = g(:,1)+X*l(1);
g(:,2) = g(:,2)+Y*l(2);
g(:,2:3) = ((n+1)*l(2)/(2*pi))*[sin(2*pi*g(:,2)/((n+1)*l(2))), cos(2*pi*g(:,2)/((n+1)*l(2)))]; % amplituda je (n+1)*l(2), delim ji pi -> prumer, ten jeste delim dvema, protoze gon. fce jdou od -1 do 1 --> v absolutni hodnote 2

% plot3(g(:,1),g(:,2),g(:,3),'r.')
l = [l(1)*(m+1) l(2)*(n+1)/pi l(2)*(n+1)/pi];
end
