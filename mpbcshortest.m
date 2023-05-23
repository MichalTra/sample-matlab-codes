function [v,d,ind] = mpbcshortest(atom1,atom2,l,overlap,once,vmat,uptriang)
% find shortest vector and distance from atom1 [x y z] to atom2 [x' y' z'] or to its image
% in the periodic boundary conditions, works even for multiple atoms in atom1 and atom2

% l is lattice matrix*
% overlap - if any atom1 is identical (even in periodic sense) with some
%        atom2, there is zero distance between them, use overlap = 1 to
%        make overlap possible (default 0 -> overlaps would be neglected)
% once - search for only the closest atom2 (default), once==true returns
%        also ind - number of closest atom2 to corresponding atom1
% vmat - return v as a matrix instead of a cell (multiple vectors on one
%        row!!!) - not very handy, but POSSIBLE SPEEDUP if once==false and
%        v is unneccessary
% uptriang - calculate only diagonal and upper-triangular elements of v and d
%        (this setting makes sense only if once==false and atom1==atom2,
%        the code will check for these two conditions!). The elements under
%        diagonal will be Inf (in matrices) or undefined (in cell).
%        POUZIVAT OPATRNE! Je to speedup a pametova optimalizace pouzitelna
%        jen za velmi specifickych podminek...

% for once==true all outputs are matrices, for once==false v is a cell
%        (unless vmat==true)

% Note: if both atom1 and atom2 contain only one atom, pbcshortest.m should
% be faster

% * When overlapping two materials with similar (but not same) unit cells,
%   it may be useful to define two unit cells. Function is able to accept
%   that using l as {l1;l2}. Function cannot handle shifts in axis direction.
%   This option should be used with caution (I hope you know what you are doing)

% Using my own implementation of mat2cell. If the standard implementation get
% faster, I can use it instead (mymat2cell case,
% mymat2celluptriang will be still needed)

% last modified: 17.4.2019
% category: chemistry, math

% EXAMPLES
%{
[v,d]=mpbcshortest([1 1 1; 2 2 2],[1 1 1; 0 2 4; 9 9 9],diag([10 10 10]),0,0); % no overlap, all shortest distances
[v,d]=mpbcshortest([1 1 1; 2 2 2],[1 1 1; 0 2 4; 9 9 9],diag([10 10 10]),1,0); % overlap, all shortest distances
[v,d,ind]=mpbcshortest([1 1 1; 2 2 2],[1 1 1; 0 2 4; 9 9 9],diag([10 10 10]),1,1); % overlap, shortest distance
[v,d,ind]=mpbcshortest([1 1 1; 2 2 2],[1 1 1; 0 2 4; 9 9 9],diag([10 10 10]),0,1); % no overlap, shortest distance
%}

if ~iscell(l)
    l2 = l;
else
    l2 = l{2};
    l = l{1};
end

if nargin == 3
    overlap = false;
    once = true;
end
if nargin == 4
    once = true;
end
if nargin<6
    vmat = false; % return cell by default, not matrix (slower, but easier to process)
end
if nargin<7
    uptriang = false; % calculate whole matrix by default, not just the upper triangle
end
if uptriang==true
    if any(size(atom1)~=size(atom2)) || any(any(atom1~=atom2))
        error('Calculating only upper triangular matrix is allowed only if atom1 and atom2 are the same!');
    end
    if once == true
        error('Calculating only upper triangular matrix makes sense only for once == false, i.e., when we are searching for all shortest contacts');
    end
end

if isempty(atom1) || isempty(atom2)
    disp('Warning: mpbcshortest: at least one matrix of coordinates (atom1, atom2) is empty!');
    d = [];
    if once
        v = [];
    else
        v = {};
    end
    ind = [];
    return
end

% move all atoms inside unit cell
atom1 = atom1 - floor(atom1/l)*l;
atom2 = atom2 - floor(atom2/l2)*l2;

% meshgrid slow -> writing meshgrid(-1:1, -1:1, -1:1) explicitly
X = [-1    -1    -1; -1     0    -1; -1     1    -1
      0    -1    -1;  0     0    -1;  0     1    -1
      1    -1    -1;  1     0    -1;  1     1    -1
     -1    -1     0; -1     0     0; -1     1     0
      0    -1     0;  0     0     0;  0     1     0
      1    -1     0;  1     0     0;  1     1     0
     -1    -1     1; -1     0     1; -1     1     1
      0    -1     1;  0     0     1;  0     1     1
      1    -1     1;  1     0     1;  1     1     1];

a2 = atom2';
a2 = (a2(:))'; % atoms 2 are in row now [x1 y1 z1 x2 y2 z2]

a2 = repmat(X*l2,1,size(atom2,1)) + repmat(a2,size(X,1),1);

if once
    v = Inf(size(atom1,1),3); % Inf preallocated - do not change, it is signalling empty data
    d = Inf(size(atom1,1),1); % Inf preallocated - do not change, it is signalling empty data
    ind = zeros(size(atom1,1),1);
else
    v = Inf(size(atom1,1),3*size(atom2,1)); % Inf preallocated - do not change, it is signalling empty data
    d = Inf(size(atom1,1),size(atom2,1)); % Inf preallocated - do not change, it is signalling empty data
end

if ~uptriang
    for i=1:size(atom2,1)
        D = distmat(atom1,a2(:,(i-1)*3+1:i*3));
        if ~overlap
            D(D==0) = Inf;
        end
        [pos,hod] = locateextreme(@min,D,1);
        if once
            zapsat = d>hod; % index for writing
            d(zapsat) = hod(zapsat);
            v(zapsat,:) = a2(pos(zapsat,2),(i-1)*3+1:i*3)-atom1(zapsat,:);
            ind(zapsat) = i;
        else
            d(:,i) = hod;
            v(:,(i-1)*3+1:i*3) = a2(pos(:,2),(i-1)*3+1:i*3)-atom1;
        end
    end
else % upper triangular matrix only (ONCE MUST BE FALSE)
    for i=1:size(atom2,1)
        D = distmat(atom1(1:i,:),a2(:,(i-1)*3+1:i*3));
        if ~overlap
            D(D==0) = Inf;
        end
        [pos,hod] = locateextreme(@min,D,1);
        d(1:i,i) = hod;
        v(1:i,(i-1)*3+1:i*3) = a2(pos(:,2),(i-1)*3+1:i*3)-atom1(1:i,:);
    end
end

if ~once && ~vmat
    % v = mat2cell(v,ones(size(v,1),1),3*ones(size(v,2)/3,1));
    if ~uptriang
        v = mymat2cell(v,1,3);
    else
        v = mymat2celluptriang(v,1,3);
    end
end

end

function V = mymat2cell(v,n,m)
% standard mat2cell slow -> this replacement is not general, but faster (still quite slow)
V = cell(size(v,1)/n,size(v,2)/m);

for i = 1:size(V,1)
    for j = 1:size(V,2)
        V{i,j} = v((i-1)*n+1:i*n,(j-1)*m+1:j*m);
    end
end

end

function V = mymat2celluptriang(v,n,m)
% standard mat2cell slow -> this replacement is not general, but faster (still quite slow)
V = cell(size(v,1)/n,size(v,2)/m);

for i = 1:size(V,1)
    for j = i:size(V,2)
        V{i,j} = v((i-1)*n+1:i*n,(j-1)*m+1:j*m);
    end
end

end

function [pos,hod]=locateextreme(handle,X,one)
% provide indices of minima/maxima on a matrix row (for columns use transpose)
% Useful especially if the position of the minima/maxima is more important than the actual value

% pos ... position (2 columns)
% hod ... value

% (handle @min/@max)
% one=1 -> search only 1st occurence

% last modified: 11.6.2018 (Mat-Oct compatible)
% category: math

if one
    [hod,pos] = handle(X,[],2);
    pos = [(1:size(X,1))' pos];
    return
end
    
hod=handle(X,[],2);
ind=X-repmat(hod,1,size(X,2));
ind=ind==0;
[a,b]=find(ind);
ab=sortrows([a(:) b(:)]); % if X have only one row, returns a and b as a row, which is what I don't want -> (:)
pos=ab;
if nargout==2 && one~=1 % in case of multiple minima length(hod)<=length(pos)
    hod=X(sub2ind(size(X),pos(:,1),pos(:,2)));
end
end

function D = distmat(a,b)
% distance matrix between two matrices of points in 3D

% a and b are 3-column matrices

% last modified: 5.8.2021
% category: math

rx = a(:,1) - b(:,1)';
ry = a(:,2) - b(:,2)';
rz = a(:,3) - b(:,3)';

D = sqrt(rx.^2 + ry.^2 + rz.^2);

end
