function [a,b] = rcerovinyconvert(c,d)
% konverze rovnice roviny ve tvaru ax+by+cz+d=0 do libovolneho
% bodu a normaly roviny nebo naopak

% USAGE:
% [bod,normr] = rcerovinyconvert(rce);
% rce = rcerovinyconvert(bod,norm);

% EXAMPLES:
% [bod,normr] = rcerovinyconvert([0 0 1 0]);
% rce = rcerovinyconvert([0 0 0],[0 0 1]);

% POZNAMKY:
% Pozor na to, ze pri pouziti [bod,normr] = rcerovinyconvert(rce) funkce
% vraci nejaky nahodny bod -> opakovane pouziti konverze tam a zpatky
% nevraci souradnici stejneho bodu!

% last modified: 17.3.2020
% cateogry: math

if nargin == 2
    c = c(:)';
    d = d(:)';
    
    D = -dot(c,d); % D je d z analyticke rovnice roviny, zatimco d je normala roviny
    a = [d D];
elseif nargin == 1
    c = c(:)';
    
    b = c(1:3);
    % hledani libovolneho bodu
    cc = 1;
    while(true)
        xx = rand(1,2);
        p = randperm(3,2);
        m = setdiff([1 2 3],p);
        a = [0 0 0];
        a(p) = xx;
        a(m) = -(b(p)*xx' + c(4))/b(m);
        % kriterium prijeti bodu jsem musel trochu zjemnit (kvuli
        % numerickym sumum to nekdy nekonverguje do mrte -> pokud probehne
        % hodne cyklu, tak trochu povolim konvergencni kriteria)
        if b*a'+c(4)==0 || (cc>10 && abs(b*a'+c(4))<1e-20) || (cc>25 && abs(b*a'+c(4))<1e-15)
            break;
        end
        cc = cc+1;
    end
else
    error('Pro prevod normaly a bodu na analyticky tvar roviny jsou nutne dva inputy, pro prevod analytickeho tvaru na normalu a libovolny bod je potreba jeden input.');
end

end