function [fbasalt,ftonalite,dfbdT,dftdT,T]=meltfraction(Pm)

%inputs: structure Pm with parameter values
%outputs: fbasal, ftonalite - melt fractions of host and dike material
% dfbdT, dftdT - derivatives of melt fraction wrt temperature
% T - temperature vector

% Copyright (C)2018- Leif Karlstrom

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

T=700:1200;

% Pm.Tliq=1100; %tonalite liquidus (C)
% Pm.Tsol=725;
% Pm.Tliqd=1165; %basalt liquidus (C)
% Pm.Tsold=1015;

ftonalite=zeros(size(T));
fbasalt=zeros(size(T));

for i=1:length(T)
    if T(i)>=Pm.Tsol && T(i)<=Pm.Tliq
        %model tonalite as piecewise linear, based on fit to P&D2005 
        if T(i)<=800
            ftonalite(i) = (0.2/(800-Pm.Tsol))*(T(i)-Pm.Tsol);
        elseif T(i)> 800 && T(i) <= 895
            ftonalite(i) = 0.2;
        elseif T(i) > 895 && T(i) <= 950
            ftonalite(i) = 0.2+(0.9-0.2)/(950-895)*(T(i)-895);
        elseif T(i)>950
            ftonalite(i) = 0.9 + (1.0-0.9)/(Pm.Tliq-950) *(T(i)-950);
        end
    end
    
    if T(i)>=Pm.Tsold && T(i)<=Pm.Tliqd
        %model basalt as powerlaw, based on fit to P&D2005 
        fbasalt(i) = ((T(i)-Pm.Tsold)/(Pm.Tliqd-Pm.Tsold))^Pm.bd;
    end    
end

dfbdT=gradient(fbasalt,T(2)-T(1));
dftdT=gradient(ftonalite,T(2)-T(1));

dftdT(dftdT<0)=0;
dfbdT(dfbdT<0)=0;


