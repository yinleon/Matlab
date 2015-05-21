%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the Matlab code for PS #3 for Professor Olivier Pauluis'
% Atmosphere, Ocean, and Climate Dynamics Course at NYU.
% Submitted Monday March 9th 2015.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%load variables
load data.mat
%define plotting parameters
Fsize = 20;
Lwide = 4;
set(gca,'Fontsize',Fsize);
%Constants
p0 = 1000;          %mb standard ref. pressure
R = 287;            %J/mol *K gas constant
Kappa = .286;       %R/cp
g = 9.81;           %m/s acceleration due to gravity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Compute the potential temperature at each level and plot both the
%    temperature and potential temperature as function of the pressure.
%    Can you locate the tropopause in the sounding?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pres = Pressuremb1;                     %Pressure in millibars
TempK = TemperatureC1+273.15;           %Kelvin
N = length(Pressuremb1);                %number of entries
pt = zeros(1,N);                        %produce empty array for calculated potential temp.
H = (R*TempK(1))/g;                     %m scale height

for i = 1:N
    pt(i) = TempK(i)*(p0/Pres(i))^Kappa;%equation for potential temperature.
end

%Plot the graph
h1=plot(pt,Pres,'Color','k','LineWidth',Lwide);
hold on
h2=plot(TempK,Pres,'Color','r','LineWidth',Lwide);
hold on
xlabel('Temperature (K)')
ylabel('Pressure (mb)')
title('Potential Pressure of Upton NY in March')
legend([h1,h2],'Potential Teperature','Observed Temperature')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Using the hydrostatic balance and ideal gas law, determine the height
%    (z) of the 500hPa and 200hPa pressure levels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1 = 500; %mb
p2 = 200; %mb
z1 = H*log(double(p0/p1));
z2 = H*log(double(p0/p2));
Z1=['The altitude of pressure at 500hPa is: ',num2str(z1),'meters'];
Z2=['The altitude of pressure at 200hPa is: ',num2str(z2),'meters'];
disp(Z1)
disp(Z2)
zMax = H*log(double(p0/Pres(N))); %test for the height at the last pressure point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Compute the Brunt-Vaisala frequency at 1000hPa, 500hPa, and 100hPa.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the indicies for pressures of interest. 
ind=find(Pres==1000);
ind1=find(Pres==500);
ind2=find(Pres==100);
% Note that the pressure and temperature arrays are parallel...

% to solve this equation the Brunt-Vaisala equation is observed in terms of
% T rather than theta, and discretized to handle the first order derivative (d(theta)/dz).
pt1000 = (g /(TempK(ind)*(p0/Pres(ind))^Kappa))*(((TempK(ind)*(p0/Pres(ind))^Kappa)-(TempK(ind)*(p0/Pres(ind-1))^Kappa))/((H*log(p0/Pres(ind))-(H*log(p0/Pres(ind-1))))))^(.5);
pt500 = (g /(TempK(ind1)*(p0/Pres(ind1))^Kappa))*(((TempK(ind1)*(p0/Pres(ind1))^Kappa)-(TempK(ind1)*(p0/Pres(ind1-1))^Kappa))/((H*log(p0/Pres(ind1))-(H*log(p0/Pres(ind1-1))))))^(.5);
pt100 = (g /(TempK(ind2)*(p0/Pres(ind2))^Kappa))*(((TempK(ind2)*(p0/Pres(ind2))^Kappa)-(TempK(ind2)*(p0/Pres(ind2-1))^Kappa))/((H*log(p0/Pres(ind2))-(H*log(p0/Pres(ind2-1))))))^(.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Was there convection that day?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ALRd = -.010; %(K/m) dry adibatic lapse rate
ALRs = -.003; %(K/m) saturated adibatic lapse rate;
TempK(1); %surface temperature

for j = 1:N
    z = H*log(double(p0/Pres(i))); %(m) altitude
    %(TempK(j)-TempK(1))/z;
    if (((TempK(j)-TempK(1))/(z)) <= ALRd) % checks instablity under dry conditons.
        disp('convection occurs at ',z,' under dry conditions.');
    end   
    if (((TempK(j)-TempK(1))/(z)) <= ALRs) % checks instability under saturated conditons.
        disp('convection occurs at ',z,' under saturated conditions.');
    end
end    
