function [x,time,theta,f] = diffusion_1d_dike_wrapper(Pm)
%(!!!!!) assumes that some of the parameters are defined outside of this wrapper 

% Copyright (C)2021 - Leif Karlstrom

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.


%% parameters that you probably don't need to change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Pm.hostrockmeltmodel == 2
Pm.b=1.7; %exponent for melt fraction temperature curve. b-->0 is "more silicic"
Pm.Tsol=1015; %host rock solidus in C
Pm.Tliq=1145; %host rock liquidus in C
elseif Pm.hostrockmeltmodel == 1
Pm.Tsol=725; %host rock solidus in C
Pm.Tliq=1060; %host rock liquidus in C
end


%boundary condition in far field away from dike. 
%RightBC=1 is insulating (Neumann condition), RightBC=2 is Dirichlet (fixed temperature at initial conditions)
Pm.RightBC=2; 

Pm.R=8.314*1e-3; %gas constant kJ/mol K
Pm.rho=2600 ;%host rock density, kg/m3
Pm.rhom=2800 ;%magma density, kg/m3
Pm.Lf=300000; %latent heat fusion J/Kg
Pm.Tsold=1015; %dike material solidus in C
Pm.Tliqd=1145; %dike material liquidus in C
Pm.bd=1.7; %exponent for basaltic dike, from MELTS calculations
Pm.cp = 1100; %specific heat capacity (assumed same between materials)

Pm.kappa = Pm.k/(Pm.cp*Pm.rho); %thermal diffusivity

%coordinate transform for increased spatial resolution near dike
%for the Currie pt study, we will want uniform grid spacing - slower, but
%better near the dike
Pm.dononuniform=0; %flag for uniform or non-uniform grid spacing
Pm.ell=Pm.L/11; %scale parameter for non-uniform grid (arctangent function) - smaller values are more uniform grid spacing
Pm.N = Pm.L/Pm.h ; %number of nodes = domain length in km/spacing


% run thermal diffusion model
tic
[x,time,theta,f] = diffusion_1d_dike(Pm);
toc

if Pm.dononuniform
    dZ=1/(Pm.N);
    dikeThickZ=atan(Pm.DikeThick/Pm.ell)/atan(Pm.L/Pm.ell);
    dikepts = 1:round(dikeThickZ/dZ);
    Z = 0:dZ:1; %distance vector in scaled coordinates
else
    dikepts = 1:round(Pm.DikeThick/Pm.h);
end    

end

