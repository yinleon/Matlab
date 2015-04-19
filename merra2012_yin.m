%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is code written by Leon Yin for an assignment in Olivier Pauluis'  %
% Atmosphere, Ocean, and Climate Dynamics course taught at Courant        %
% Institute of Mathematics.                                               %
%                                                                         %
% All referenced figures are from Marshall and Plumb's textbook of the    %
% same name.                                                              %
%                                                                         %
% This code uses reanalysis data from NASA satalite observations.         %
% Data taken from a 2012 reanalysis fron NASA's Merra system contains     %
% rich 4-dimensional(longitude,latitude,pressure,month) data fields.      %
% The data is scrubbed and used to calculate the zonally averaged         %
% (1a) temperature, (1b) wind velocity, (2) mass transport, and           %
% (3) reletive humidity for January and July 2012. The data is used       %
% to analyze convection cells, jetstreams, and atmospheric moisture.      %
%                                                                         %
% Submitted: April 20th, 2015                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% Global Variables
load MERRA2012.mat
g= 9.81;                %(m/s)  acceleration due to gravity.
omega = 1/867400;       %(s)    Earth's rotation.
a = 6378100;            %(m)    radius of earth.
r = 0.0174532925;       %       radians per degree.
Epsilon = 0.62198;      %       molecular weights saturated:dry air.
lLat = length(lat);     %       dimension of latitude.
lP = length(plev);      %       dimension of pressure.
Jan = 1;
July = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: The format of data is as follows:                                 %
%                                                                         %
%       PS(lat, month) surface pressure                                   %
%       U(lat, pressure, month) mean zonal wind                           %
%       V(lat, pressure, month) mean meridional wind                      %
%       Q(lat, pressure, month) mean specific humidity                    %
%       T(lat, pressure, month) mean temperature                          %
%       Z(lat, pressure, month) mean geopotential height                  %
%                                                                         %
%       plev(pressure) pressure at each level (in hPa)                    %
%       lat(lat) latitude (in degree)                                     %
%                                                                         %
%       The dimension of the data is as follows:                          %
%       144 latitude                                                      %
%       42 pressure levels                                                %
%       12 months                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data Scrubbing
% Some values in the MERRA2012.mat file are 'NaN' (not a number), and must
% be removed for the integrity of calculations.

Ujan  = zeros(lLat,lP);
Vjan  = zeros(lLat,lP);
Qjan  = zeros(lLat,lP);
Tjan  = zeros(lLat,lP);
Ujuly = zeros(lLat,lP);
Vjuly = zeros(lLat,lP);
Qjuly = zeros(lLat,lP);
Tjuly = zeros(lLat,lP);

% Iterate through each dataset, and replace every instance of 'NaN'.
for L = 1:lLat
    for P = 1:lP
        % For U of Jan
        if isnan(U(L,P,Jan)) 
            Ujan(L,P) = 0.0;
        else
            Ujan(L,P) = U(L,P,Jan);
        end
         % For V of Jan
        if isnan(V(L,P,Jan)) 
            Vjan(L,P) = 0.0;
        else
            Vjan(L,P) = V(L,P,Jan);
        end
         % For T of Jan
        if isnan(T(L,P,Jan)) 
            Tjan(L,P) = 255.0;
        else
            Tjan(L,P) = T(L,P,Jan);
        end
         % For Q of Jan
        if isnan(Q(L,P,Jan)) 
            Qjan(L,P) = 0.0;
        else
            Qjan(L,P) = Q(L,P,Jan);
        end  
        % For U of July
        if isnan(U(L,P,July)) 
            Ujuly(L,P) = 0.0;
        else
            Ujuly(L,P) = U(L,P,July);
        end
        % For V of July
        if isnan(V(L,P,July))
            Vjuly(L,P) = 0.0;
        else
            Vjuly(L,P) = V(L,P,July);
        end
        % For T of July
        if isnan(T(L,P,July)) 
            Tjuly(L,P) = 255.0;
        else
            Tjuly(L,P) = T(L,P,July);
        end
         % For Q of July
        if isnan(Q(L,P,July)) 
            Qjuly(L,P) = 0.0;
        else
            Qjuly(L,P) = Q(L,P,July);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Plot the zonal mean and temperature for January and July. Compare    %
%    the location and intensity of the jet stream (i.e. the strong        %
%    westerly wind at about 200 hPa) in both seasons hemispheres.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables and Preconditions
UzJan  = transpose(Ujan);       %transpose mean zonal wind -> U(P,Lat)                                   
UzJuly = transpose(Ujuly);      %transpose mean zonal wind -> U(P,Lat)
TzJan  = transpose(Tjan);       %transpose mean zonal temp -> T(P,Lat) 
TzJuly = transpose(Tjuly);      %transpose mean zonal temp -> T(P,Lat)

% Plot the graphs for January
figure(1)
C = contourf(lat,plev,UzJan,15); 
set(gca,'ydir','rev');  %invert the y-axis to show data from ground->up.
clabel(C)
xlabel('Latitude')
ylabel('Pressure (mb)')
title('Zonal-Average Wind Velocity in January 2012')
colorbar

figure(2)
C = contourf(lat,plev,TzJan,15); 
set(gca,'ydir','rev');  %invert the y-axis to show data from ground->up.
clabel(C)
xlabel('Latitude')
ylabel('Pressure (mb)')
title('Zonal-Average Temperature in January 2012')
colorbar

% Plot the graphs for July
figure(3)
C = contourf(lat,plev,UzJuly,15);
set(gca,'ydir','rev');  %invert the y-axis to show data from ground->up.
clabel(C)
xlabel('Latitude')
ylabel('Pressure (mb)')
title('Zonal-Average Wind Velocity in July 2012')
colorbar

figure(4)
C = contourf(lat,plev,TzJuly,15); 
set(gca,'ydir','rev');  %invert the y-axis to show data from ground->up.
clabel(C)
xlabel('Latitude')
ylabel('Pressure (mb)')
title('Zonal-Average Temperature in January 2012')
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Compute the streamfunciton for January and February associated with  %
%    the meridional circulation (figure 5.21 from the textbook).          %
%                                                                         %
%    The Streamline function is as follows:                               %
%           Ψ(lat,P) = 2*pi*a*cos(lat)*(1/g)*integral of (Vm dp)          %
%    Where:                                                               %
%           P  is pressure,                                               %
%           a  is the radius of earth,                                    %
%           Vm is meridional velocity.                                    %
%                                                                         %
%    The integral can be solved implicitly with trapezoidal rule:         %
%           (P2-P1)*(V(P1) + V(P2))/2                                     %
%                                                                         %
%    NOTE: cos(arg) takes an argument in radians.                         %
%                                                                         %
%    Identify the Hadley and Ferrel cells in the figures and indicate     %
%    the ascent and subsidence regions. Discuss the intensity of these    %
%    circulations and how they vary between the winter and summer.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables and Preconditions
streamJan  = zeros(lLat,lP);
streamJuly = zeros(lLat,lP);

% Calculations
for L = 1:lLat-1
    % Calculate the Reiman Sum for each latitude point.
    reimanSumJan = 0;
    reimanSumJuly = 0;
    for P = 1:lP-1
        % Increment ReimanSumJan per pressure level.
        reimanSumJan = reimanSumJan...
            +(((Vjan(L,P)+Vjan(L,P+1))/2)*(plev(P+1)-plev(P)));
        % Calculate Stream Function for January.
        streamJan(L,P) =...
            (2 * pi * a * cos(lat(L)*r) * (1/g)) * reimanSumJan;
        
        % Increment the ReimanSumJuly per pressure level.
        reimanSumJuly = reimanSumJuly...
            +(((Vjuly(L,P)+Vjuly(L,P+1))/2)*(plev(P+1)-plev(P)));
        % Calculate Stream Function for July.
        streamJuly(L,P) = ...
            (2 * pi * a * cos(lat(L)*r) * (1/g)) * reimanSumJuly;
    end
end

% Plot the graph for January
figure(5)
C = contourf(lat,plev,transpose(streamJan),20);
set(gca,'ydir','rev');  %invert the y-axis to show data from ground->up.
clabel(C)
xlabel('Latitude')
ylabel('Pressure (mb)')
title('Zonal-Average Streamfunction January 2012')
colorbar

% Plot the graph for July
figure(6)
C = contourf(lat,plev,transpose(streamJuly),20);
set(gca,'ydir','rev');  %invert the y-axis to show data from ground->up.
clabel(C)
xlabel('Latitude')
ylabel('Pressure (mb)')
title('Zonal-Average Streamfunction July 2012')
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Compute the relative humidity based on expression (4-25) for January %
%    and July. Where are the dryer and moister regions of the atmosphere? %
%    Can you relate these to the circulation computed in (2)?             %
%                                                                         %
%    Saturated Humidity can be calculated as follows:                     %
%           humSat = Epsilon*Es/p              (4-23)                     %
%    Where:                                                               %
%           Epsilon is the ratio of saturated:dry air molecular weight.   %
%           Es      is saturated vapor pressure                           %
%                                                                         %
%    Using this calc, and specific humidity Q, we can get relative        %
%    humidity using the following equation:                               %
%           humRel = Q/humSat * 100            (4-25)                     %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables and Preconditions
humRelJuly = zeros(lLat,lP);
humRelJan  = zeros(lLat,lP);

% Calculations
for L = 1:lLat
    for P = 1:lP
        % Calculate saturated vapor pressure for July.
        eSJuly = 6.11*exp(0.067*(Tjuly(L,P)-273.15)); 
        % Calculate saturated humidity for July.
        humSat = Epsilon*eSJuly/plev(P);
        % Calculate relative humidity of July.
        humRelJuly(L,P) = (Qjuly(L,P)/humSat)*100;
        
        % Calculate saturated vapor pressure for Jan.
        eSJan  = 6.11*exp(0.067*(Tjan(L,P)-273.15));
        % Calculate saturated vapor pressure for Jan.
        humSat = Epsilon*eSJan/plev(P);    
        % Calculate relative humidity of Jan.
        humRelJan(L,P) = (Qjan(L,P)/humSat)*100;
    end
end

% Plot the graph for January
figure(7)
C = contourf(lat,plev,transpose(humRelJan),15);
set(gca,'ydir','rev');  %invert the y-axis to show data from ground->up.
clabel(C)
xlabel('Latitude')
ylabel('Pressure (mb)')
title('Zonal-Average Relative Humidity (%) January 2012')
colorbar

% Plot the graph for July
figure(8)
C = contourf(lat,plev,transpose(humRelJuly),15);
set(gca,'ydir','rev');  %invert the y-axis to show data from ground->up.
clabel(C)
xlabel('Latitude')
ylabel('Pressure (mb)')
title('Zonal-Average Relative Humidity (%) July 2012')
colorbar
