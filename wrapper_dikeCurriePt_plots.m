%wrapper for 1D thermal model in Biasi and Karlstrom (2021), "Timescales of 
%magma transport in the Columbia River flood basalts, determined by paleomagnetic 
%data" in revision at Earth and Planetary Science Letters

%Use this code to generate figures with either 2D parameter sweep or 
%multiple curves varying one parameter. Will generate results equivalent to
%figures in Supplement.

% Copyright (C)2021 - Leif Karlstrom

% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BASIC SWITCH: either plot a few curves or a 2D parameter sweep
plottype = 1; %plot a few curves of Curie temperature through time varying a single parameter
%plottype = 2; %plot a 2D parameter sweep with max Curie distance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch plottype
    case 1
        %define a vector for whatever parameter you want too vary. You'll
        %need to modify the input data structure accordingly if you change
        %ParameterVec = [10]; %this is thermal conductivity in W/mC
        %ParameterVec = [1 2]; %this is host rock composition
        %ParameterVec = [25 50 75 100]; %this is background temperature in C
        %ParameterVec = [0.5 5 10]; %this is dike half-thickness in m
        ParameterVec = [10]; %this is flow unsteadiness Tw in years
        %ParameterVec = [10]; %this is dike active lifetime Tc in
        %years
        
        
        for ii = 1:length(ParameterVec)
            disp(ii)
            
            clear x time theta f 

            Pm.year = 3600*24*365;  %one year in seconds            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %things you might want to vary
            Pm.hostrockmeltmodel = 1; %ParameterVec(ii);%1; %1=tonalite, 2=basalt
            Pm.k = 3; %ParameterVec(ii); %2; %thermal conductivity (assumed same between materials)
            Pm.twidth = ParameterVec(ii); %ParameterVec(ii); 0.01; %vary this parameter to change the rate at which dike temperatures decrease - larger values are less sharp step
            Pm.tcenter = 1; %ParameterVec(ii); %10*3600*24/Pm.year; %scale for dike shut-off time (units of years)
            Pm.Tbackground = 40; %ParameterVec(ii); %40; %background temperature in C
            Pm.DikeThick = 5; %ParameterVec(ii);%5; %dike half-thickness in m
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %things that control resolution and domain size (good to test
            %sensitivity of results to changing these)
            Pm.h = 0.1;%Pm.L/100; %grid spacing, m
            Pm.L = 30; %total length of domain in meters - can be relatively short since we're focused on very high temperature constraint
            Pm.TotT = 10;  %total time of simulation in years - can also be relatively short due to high temp
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %run the thermal diffusion model (including wrapper)
            [x,time,theta,f] = diffusion_1d_dike_wrapper(Pm);
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %post-processing and plotting
            
            %find the 580 deg isotherm
            Currie = 580;
            [~,Cval]=min(abs(max(theta,[],2)-Currie));

            if(max(theta(Cval-1,:))>Currie && max(theta(Cval,:))<Currie)
                Xcmin(ii)=x(Cval-1) - Pm.DikeThick; Xcmax(ii)=x(Cval)- Pm.DikeThick; 
            else
                Xcmin(ii)=x(Cval) - Pm.DikeThick; Xcmax(ii)=x(Cval+1)- Pm.DikeThick;
            end
            
            XaboveC=zeros(size(time));
            %then as a timeseries
            for jj=1:length(time)
                TaboveC= sum(theta(:,jj)>Currie);
                if TaboveC>0
                XaboveC(jj) = x(TaboveC);
                end
            end
            
            figure(3)
            plot(time/Pm.year,XaboveC,'o-');hold on
            xlabel('Time in years'); 
            ylabel('Distance from dike center above 580 C')
            labl3{ii,:} = ['Parametervalue = ' num2str(ParameterVec(ii))];
            
        end
            
            %plot up a figure illustrating variation in Currie point temp
            %as fct of parameter
            figure(1)
            for i=1:length(ParameterVec)
                plot([ParameterVec(i) ParameterVec(i)],[Xcmin(i) Xcmax(i)],'b-','LineWidth',5); 
                hold on; 
            end
            xlabel('Parameter value');
            ylabel('Max distance above 580 C from dike wall (based on grid res.)')
            
            
            %plot up the temperature as a function of distance from dike
            %center, at several times to illustrate thermal evolution
            
            Tvec=[0 .001 .01 .1 .2 .5 1 1.5 2 3 4 5 6 7 8 9 10];%linspace(0,Pm.TotT,int);

            
            figure(2)
            for i=1:length(Tvec)
                [~,tval]=min(abs(time-Tvec(i)*Pm.year));
            plot(x,theta(:,tval),'-')
            hold on
            lablv{i,:} = ['T = ' num2str(Tvec(i)) ' years'];
            end
            xlabel('distance from dike center in meters')
            ylabel('Temperature in C')
            xlim([0 30])
            legend(lablv)        

            figure(3)
            legend(labl3)
            

    case 2
        %assume that we want to explore Tcenter and Twidth to make a
        %parameter plot. Will pick "typical" values for k and dike
        %thickness
        
                Pm.year = 3600*24*365;  %one year in seconds
                
                lenTc = 30; %number of steps in Tc
                lenTw = 35; %number of steps in Tw
                
                Tc = logspace(0.3,3.27,lenTc)*3600*24/Pm.year;%scale for dike shut-off time (years)
                Tw = logspace(-2,1,lenTw);%scale for unsteadiness in dike boundary temperature (years)
            	
                
                for ii = 1:length(Tc)
                    for jj = 1:length(Tw)
                        
                    clear x time theta f 
    
                        
                        %time (in years) at which tanh reaches 0.01 of initial value. This is total duration of the active dike.
                        Tf(ii,jj) = atanh(1-0.01*2)*Tw(jj)*0.1 + Tc(ii); 
                        
                        disp(ii)

                        %%%%%%%%%%%%%%%%%%%%%%%%
                        %things you might want to vary
                        Pm.hostrockmeltmodel = 1; %1=tonalite, 2=basalt
                        Pm.k = 3;%ParameterVec(ii); %2; %thermal conductivity (assumed same between materials)
                        Pm.Tbackground = 40; %background temperature in C
                        Pm.DikeThick = 5; %dike half-thickness in m

                        Pm.twidth = Tw(jj); %vary this parameter to change the rate at which dike temperatures decrease - larger values are less sharp step
                        Pm.tcenter = Tc(ii); %scale for dike shut-off time (years)

                        %%%%%%%%%%%%%%%%%%%%%%%%
                        %things that control resolution and domain size (good to test
                        %sensitivity of results to changing these)
                        Pm.h = 0.2;%Pm.L/100; %grid spacing, m
                        Pm.L = 30; %total length of domain in meters - can be relatively short since we're focused on very high temperature constraint
                        Pm.TotT = 10; %min(6*Tf(ii,jj),75);%25;  %total time of simulation in years - can also be relatively short due to high temp

                        %%%%%%%%%%%%%%%%%%%%%%%%
                        %run the thermal diffusion model (including wrapper)
                        [x,time,theta,f] = diffusion_1d_dike_wrapper(Pm);

                        %%%%%%%%%%%%%%%%%%%%%%%%
                        %post-processing and plotting

                        %find the 580 deg isotherm
                        Currie = 580;
                        [~,Cval]=min(abs(max(theta,[],2)-Currie));

                        if(max(theta(Cval-1,:))>Currie && max(theta(Cval,:))<Currie)
                            Xcmin=x(Cval-1) - Pm.DikeThick; Xcmax=x(Cval)- Pm.DikeThick; 
                        else
                            Xcmin=x(Cval) - Pm.DikeThick; Xcmax=x(Cval+1)- Pm.DikeThick;
                        end

                        %choose the minimum point as conservative estimate of Tc
                        %distance for MGT
                        MaxTc(ii,jj) = max(0,Xcmin);
                    end
                end
                
                
                %% 
                
                %make a contour plot of the maximum reset distance as a
                %function of Tw and Tw, overlay with contours of the total
                %active lifetime Tf
                
                 figure(1)   
                 %Cont=[0 0.18 0.2 0.23 0.26 0.28 0.6 0.65 0.66 1 2 3 4 5 6 7 8 9 10]; %defince contours of MaxTc (in m) you want to see
                 Cont=linspace(0,10);
                 contourf(Tw,Tc,MaxTc,Cont)
                 set(gca,'xscale','log')
                 set(gca,'yscale','log')
                 a = linspace(0.25,1);
                 b = linspace(0.25,0);
                 c = linspace(1,0);
                 map = [a',b',c'];
                 colormap(map)
                 colorbar
                 xlabel('Tw (years)')
                 ylabel('Tc (years)')
                 
                 %overlay plot of total active timescale Tf
                 TfVec=[0.03 0.1 0.25 0.5 1 2 3 4 5 6 7 8 9 10 20 30]; %vector in years
                 %TfVec=logspace(-2,1); %vector in years
                 hold on
                 contour(Tw,Tc,Tf,TfVec,'w--','ShowText','on')
                 caxis([0 10])
                 title('Colored contours = distance to 580 C in m')
                 hold off
                 
                 %%
                 %make a plot of the "magma flux" curves associated with
                 %the range of Tc and Tw used
                
                 %compare to Fig 5 and Fig12A of Karlstrom et al 2019
                 numcurves=9;
                 %TcV=logspace(log10(min(Tc)),log10(max(Tc)),3);
                 %TwV=logspace(log10(min(Tw)),log10(max(Tw)),3);
                 TcV=[1 1 1 1 1];
                 TwV=[0.01 0.1 1 2 10];
                 
                 
                 
                 figure(2)
                 num=1;
                 for j=1:length(TcV)
                     for i=1:length(TcV)
                        %define the normalized boundary temperature
                        Magma = 0.5*(1-tanh((time/Pm.year-TcV(j))/TwV(i)));
                 
                        plot((time/Pm.year),Magma); hold on
                        lablT{num,:} = ['Tc=' num2str(TcV(j)) ' yr, Tw=' num2str(TwV(i)) ' yr'];
                        num=num+1;
                     end
                 end
                 xlim([0 10])
                 legend(lablT)
                 hold off
end     


