function x = prunikprimkyaroviny(bodp,normp,bodr,normr,eos)
% najde prunik primky a roviny (hleda bod, pokud by lezela primka v rovine,
% kod havaruje)

% bodp a normp je bod a normalovy vektor definujici primku
% bodr a normr je bod a vektor definujici rovinu
% eos (error on subset) - parametr urcujici, zda ma funkce spadnout, pokud
%     je primka podmnozinou roviny (eos=true, default) nebo jestli ma
%     vratit NaN (eos=false)

% last modified: 11.3.2019
% category: math

if nargin<5
    eos = true;
end

% pokud je primka rovnobezna s rovinou
if dot(normp,normr)==0
    [~,d] = vzdalenostboduodroviny(bodp,bodr,normr);
    if d==0
        if eos
            error('Primka lezi v rovine. Prunikem je tudiz cela primka!');
        else
            x = NaN;
            return;
        end
    else
        x = [];
        return
    end
end

d = -dot(bodr,normr); % rovina ma rovnici normr(1)*x+normr(2)*y+normr(3)*z+d = 0

k = -[bodp 1]*[normr(:);d]/dot(normp,normr);

x = bodp + k*normp;

end
