function [q1,med,q3] = boxplot(data_Y,x,w,type,color)
% boxplot function

% INPUTS
% data_Y ... data
% x      ... shift in x axis (default 0)
% w      ... width
% type   ... type (0-default, 1-include outliers)
% color  ... color (default is black & white = 'kw')

% last modified: 16.8.2016
% plotting

% Example:
% boxplot(rand(100,1));

if nargin<2
    x = 0;
end
if nargin<3
    w = 0.3;
end
if nargin<4
    type = 0;
end
if nargin<5
    color = 'kw'; % 1st color border, 2nd color fill
end

y = sort(data_Y(:));
if length(y)<4
    error('boxplot.m: error: not enough input data - minimum is four values')
end

col = color;

if mod(length(y),2)==1
    med = y(ceil(length(y)/2));
else
    med = mean([y(length(y)/2) y(length(y)/2+1)]);
end

% there are more ways how the quartils may be calculated -> using method 3 from https://en.wikipedia.org/wiki/Quartile
n = floor(length(y)/4);
switch(mod(length(y),4))
    case 0
        q1 = mean([y(length(y)/4) y(length(y)/4+1)]);
        q3 = mean([y(3*length(y)/4) y(3*length(y)/4+1)]);
    case 1
        q1 = 0.25*y(n) + 0.75*y(n+1);
        q3 = 0.75*y(3*n+1) + 0.25*y(3*n+2);
    case 2
        q1 = y(ceil(length(y)/4));
        q3 = y(ceil(3*length(y)/4));
    case 3
        q1 = 0.75*y(n+1)+0.25*y(n+2);
        q3 = 0.25*y(3*n+2)+0.75*y(3*n+3);
end

rectangle('Position',[x-w/2,q1,w,med-q1],'Edgecolor',color(1),'Facecolor',color(2));
hold on
rectangle('Position',[x-w/2,med,w,q3-med],'Edgecolor',color(1),'Facecolor',color(2));
if type == 0
    plot([x x],[q3 max(y)],[col(1),'-']);
    plot([x-w/2 x+w/2],[max(y) max(y)],[col(1),'-']); % horizontal line
    plot([x x],[min(y) q1],[col(1),'-']);
    plot([x-w/2 x+w/2],[min(y) min(y)],[col(1),'-']); % horizontal line
else
    IQR = q3-q1; % Interquartile range
    if any(y>q3+1.5*IQR)
        plot([x x],[q3 q3+1.5*IQR],[col(1),'-']);
        plot([x-w/2 x+w/2],[q3+1.5*IQR q3+1.5*IQR],[col(1),'-']); % horizontal line
        bpplot(x,y(y>q3+1.5*IQR),[col(1),'x']);
        if col(1)=='w'
            disp('boxplot.m: Warning: white color border may cause that outliers won''t be visible!');
        end
    else % no outliers
        plot([x x],[q3 max(y)],[col(1),'-']);
        plot([x-w/2 x+w/2],[max(y) max(y)],[col(1),'-']); % horizontal line
    end
    if any(y<q1-1.5*IQR)
        plot([x x],[q1-1.5*IQR q1],[col(1),'-']);
        plot([x-w/2 x+w/2],[q1-1.5*IQR q1-1.5*IQR],[col(1),'-']); % horizontal line
        bpplot(x,y(y<q1-1.5*IQR),[col(1),'x']);
        if col(1)=='w'
            disp('boxplot.m: Warning: white color border may cause that outliers won''t be visible!');
        end
    else % no outliers
        plot([x x],[min(y) q1],[col(1),'-']);
        plot([x-w/2 x+w/2],[min(y) min(y)],[col(1),'-']); % horizontal line
    end
end    
end

function bpplot(x,y,format)
% plotting function

plot(x*ones(length(y),1),y,format)

end
