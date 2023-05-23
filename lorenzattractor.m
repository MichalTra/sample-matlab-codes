function [z,F]=lorenzattractor(r,s,b,start,timestep,time)
% Lorenz attractor - theory of chaos
% without arguments uses r, sigma, b equal to 28, 10, 8/3 (respectively)
% start point [1,0,0], timestep 0.001 s and time 25 s

% it seems that with timestep 0.01 s the final plot is different, so it may
% be numerically unstable

% category: math, physics

% Example:
% lorenzattractor;

if nargin==0
   r=28; s=10; b=8/3; start=[1 0 0]; time=25; timestep=0.001;
elseif nargin==3
   start=[1 0 0];
   time=25;
   timestep=0.001;
end

n=fix(time/timestep);

z=zeros(n+1,3);
z(1,:)=start;

j=1;
if n>1000   % prelocation
    F(1000)=struct('cdata',[],'colormap',[]);
else
    F(n)=struct('cdata',[],'colormap',[]);
end

for i=1:n
   z(i+1,1)=z(i,1)+s*(z(i,2)-z(i,1))*timestep;
   z(i+1,2)=z(i,2)+(z(i,1)*(r-z(i,3))-z(i,2))*timestep;
   z(i+1,3)=z(i,3)+(z(i,1)*z(i,2)-b*z(i,3))*timestep;
   if n<1000 || i==n || mod(i,n/1000)==0 % no more than 1000 points (final point added)
        plot3(z(1:i+1,1),z(1:i+1,2),z(1:i+1,3));
        F(j)=getframe;
        j=j+1;
   end
end

% movie(F,1,10);  % may be redrawn using this code
end
