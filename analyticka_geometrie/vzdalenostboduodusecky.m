function [v,d] = vzdalenostboduodusecky(bod,u1,u2)
% urci nejkratsi vektor od usecky k bodu(m) a jeho velikost, usecka je
% specifikovana body u1 a u2

% predpoklada data v radcich, tzn. bod(nx3), u1(1x3), u2(1x3)

% last modified: 14.3.2019
% category: math

% need access to vzdalenostboduodprimky.m

n = size(bod,1);

u = u2-u1; % vektor usecky
r1 = [u -u*u1']; % jedna rovina urcujici poloprostor
r2 = [-u u*u2']; % druha rovina urcujici poloprostor

v = zeros(n,3);
d = zeros(n,1);

ind =  r1*[bod'; ones(1,n)]>0 & r2*[bod'; ones(1,n)]>0; % pokud jsou body v poloprostorech danymi useckou, muzu vzdalenost bodu spocitat pomoci vzdalenosti od primky
if any(ind)
    [v(ind,:),d(ind)] = vzdalenostboduodprimky(bod(ind,:),u1,u);
end

% ted doplnim to, co se neda spocitat pomoci vzdalenosti od primky
v1 = [bod(~ind,1)-u1(1), bod(~ind,2)-u1(2), bod(~ind,3)-u1(3)];
v2 = [bod(~ind,1)-u2(1), bod(~ind,2)-u2(2), bod(~ind,3)-u2(3)];

d1 = sqrt(sum(v1.^2,2));
d2 = sqrt(sum(v2.^2,2));

vdop = v2;
vdop(d1<d2,:) = v1(d1<d2,:);
v(~ind,:) = vdop;
d(~ind) = min(d1,d2);
