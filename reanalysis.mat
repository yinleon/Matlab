%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is code written by Leon Yin for an assignment in Olivier Pauluis' %
% Atmosphere, Ocean, and Climate Dynamics course taught at Courant       %
% Institute of Mathematics.                                              %
% All referenced figures are from Marshall and Plumb's textbook of the   %
% same name.                                                             %
%                                                                        %
% Submitted: March 27th, 2015                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% Global Variables
load ncep_2008_clean.mat
Pause = 4;                  %(s) time between graphs of each question
p0 = 1000.;                 %(mb) reference pressure
kappa = .286;               %R
R = 287;                    %(J/mol*K) gas constant
g = 9.81;                   %(m/s) accelerationdue to gravity
T1= 250;                    %(K) averagetemperature of the atm        
H = kappa*T1/g;             %(km) scale height

% NOTE: the format for Temp is T(longitude,latitude,pressure,month)
% NOTE: the same format applies to the Geospatial Height, z_global(...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 1a) Plot a sounding of temperature and potential temperature in        %
%     New York in June. Compare to Figure 4.9 in textbook.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables and Preconditions
NYClat = 40.;               %north
NYClon = 75.;               %west
indLat=find(lat==NYClat);   %find the index for NYC latitude
indLon=find(lon==NYClon);   %find the index for NYC longitude
month = 6;                  %6th month IE June
lP = length(p);             %number of indices in the pressure array
pT = zeros(lP,1);           %creates an array for the potenetial temperature
tNY = zeros(lP,1);          %creates an array for the temperature
z = zeros(lP,1);            %creates an array for the altitude            

% Calculations
for i = 1:lP
    tNY(i) = T(indLon,indLat,i,month);      %(K) temp
    pT(i)= tNY(i)*(p0/p(i))^kappa;          %(K) potential temp
    z(i) = H*log(p0/p(i));                  %(Km) altitude
end

% Plotting parameters
Fsize = 12;
Lwide = 4;
set(gca,'Fontsize',Fsize);

% Plot the graph
h1=plot(pT,p,'Color','k','LineWidth',Lwide);
hold on
h2=plot(tNY,p,'Color','r','LineWidth',Lwide);
hold on
set(gca,'ydir','rev')     %invert the y-axis to show data from ground->up
xlabel('Temperature (K)')
ylabel('Pressure (mb)')
title('Potential Pressure NYC June 2008')
legend([h1,h2],'Potential Teperature','Observed Temperature')
pause(Pause);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1b) Is the atmosphere stable to dry adiabatic processes?               %
%     Note: if (d(pT)/d(z) < 0) IE if potential temperature decreases    %
%     with height, then the parcel is UNSTABLE and convection occurs.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations
for k=2:lP 
    if((pT(i)-pT(i-1))/(z(i)-z(i-1)) < 0)
        disp('there is convection at ',z(i),'km under dry conditions.');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Plot zonally averaged annual mean temperature.                      %
%    Compare to figure 5.7 in textbook.                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% Calculations
Tz = (mean(T,1));                       %find mean of longitude (zonal format)
TzT = transpose(squeeze(mean(Tz,4)));   %find mean per month, squeeze longitude, 
                                        %transpose data -> Temp(Pressure,Latitude)
% Plot the graph
C = contourf(lat,p,TzT,15);
set(gca,'ydir','rev');  %invert the y-axis to show data from ground->up.
clabel(C)
xlabel('Latitude')
ylabel('Pressure (mb)')
title('Zonally averaged annual mean temperature')
pause(Pause);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Plot zonally averaged annual mean potential temperature.            %
%    Compare to figure 5.8 in textbook.                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% Variables and preconditions
lLat = length(lat);
pT2 = zeros(lP,lLat);   %create an array for mean potential temperature

% Calculations
for j =1:lLat                                    %per latitude
    for i = 1:lP                                 %per pressure
        pT2(i,j)= TzT(i,j)*(p0/p(i))^kappa;      %potential temperature
    end
end

% Plot the graph
C = contourf(lat,p,pT2,25);
set(gca,'ydir','rev');  %invert the y-axis to show data from ground->up.
colormap summer;
clabel(C)
xlabel('Latitude')
ylabel('Pressure (mb)')
title('Zonally averaged annual mean potential temperature')
pause(Pause);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Plot the 500mbar geopotential height for January 2008.              %
%    Compare to figure 5.12 in textbook                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% Variables and preconditions
iP = find(p==500);     %find the index where pressure is 500 mb
month = 1;             %1st month IE January
lLon = length(lon);
geoZ=zeros(lLat,lLon); %create an array for the geopotential height at each lat and long

% Calculations
for i = 1:lLon                                  %per longitude
    for j = 1:lLat                              %per latitude
        geoZ(j,i)=z_global(i,j,iP,month);       %geopotential height
    end
end

%Plot the graph
C = contourf(lat,lon,transpose(geoZ),20);
%clabel(C)         %contour labels (optional).
xlabel('Latitude')
ylabel('Longitude')
title('Mean Geopotential height January 2008')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

