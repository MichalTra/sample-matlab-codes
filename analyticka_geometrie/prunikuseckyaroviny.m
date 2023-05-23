function x = prunikuseckyaroviny(bod1u,bod2u,bodr,normr,eos)
% prunik usecky a roviny (bod1u,bod2u,bodr,normr); hleda bod, v pripade ze
% usecka lezi v rovine, kod havaruje

% bod1u je 1. bod usecky
% bod2u je 2. bod usecky
% bodr je bod a normr normalovy vektor spolecne definujici rovinu
% eos (error on subset) - parametr urcujici, zda ma funkce spadnout, pokud
%     je usecka podmnozinou roviny (eos=true, default) nebo jestli ma
%     vratit NaN (eos=false)

% last modified: 11.3.2019
% category: math

% need access to prunikprimkyaroviny.m

% EXAMPLES:
% prunikuseckyaroviny([0 0 0],[1 1 1],[0 0.5 0],[1 2 3])

if nargin<5
    eos = true;
end

normp = bod2u - bod1u; % normala pomocne primky

x = prunikprimkyaroviny(bod1u,normp,bodr,normr,false);
if isnan(x)
    if eos
        error('Prunikprimkyausecky.m: Error, volana pomocna funkce prunikprimkyaroviny havarovala. Pravdepodobnym duvodem je, ze zadana usecka lezi v dane rovine!');
    else
        x = NaN;
        return;
    end
end

if isempty(x) % nema-li prunik primka, jiste ho nema ani usecka
    return;
end

v1 = x-bod1u;
v2 = x-bod2u;

% pokud x lezi na usecce, mely by vektory v1 a v2 byt nesouhlasne
% rovnobezne, a tudiz soucin (po slozkach) jejich znamenek by mel byt -1
% (popripade 0 v pripade, ze bod x je totozny s krajnim bodem usecky)
if ~all(sign(v1).*sign(v2)<=0)
    x = [];
end

end
