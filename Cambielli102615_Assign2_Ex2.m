%% Spacecraft Guidance and Navigation - Periodic Orbit (2023/2024)
% Assignment:     2
% Exercise:       2 - Batch filters
% Author:         Alessandro Cambielli

% ID number:      10619322
%                 102615

%% REQUEST 1

% Clear memory workspace and set path
clearvars; close all; clc 

% Load SPICE kernels: 
cspice_furnsh('assignment02.tm');
format long g

t0 = '2010-08-12-05:30:00.000 UTC'; 
tf = '2010-08-12-11:00:00.000 UTC'; 
et0 = cspice_str2et(t0);  % Ephemeris time [s]
etf = cspice_str2et(tf);  % Ephemeris time [s]
etsep = cspice_str2et('2010-08-12-05:27:39.114 UTC');

mu = cspice_bodvrd('Earth','GM',1);  % Earth gravitational parameter [km^3/s^2]
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

r0_sep = [4622.232026629; 5399.3369588058; -0.0212138165769957];
v0_sep = [0.812221125483763; -0.721512914578826; 7.42665302729053];
x0_sep = [r0_sep;v0_sep];

% Define stations name 
station1Name = 'KOUROU';
station2Name = 'SVALBARD';

[~, xx] = ode113(@(t,x) kepler(t,x,mu), [etsep et0], x0_sep, options);
x0 = xx(end,:);

% Equally spaced vector, with one minute time-step
tvec = et0:60:etf;

[~, mangostate] = ode113(@(t,x) kepler(t,x,mu), tvec, x0, options);

[azimuth1,elevation1,~,~] = antenna_pointing(station1Name, tvec, mangostate);
[azimuth2,elevation2,~,~] = antenna_pointing(station2Name, tvec, mangostate);

minelevation1 = 10; % deg
minelevation2 = 5; % deg

index_vis1 = elevation1 >= minelevation1; % KOUROU visibility window
index_vis2 = elevation2 >= minelevation2; % SVALBARD visibility window

tstart1 = [tvec(197), tvec(298)];
tend1 = [tvec(205), tvec(301)];

tstart2 = [tvec(15), tvec(117), tvec(219), tvec(321)];
tend2 = [tvec(25), tvec(125), tvec(225), tvec(326)];

fprintf('   KOUROU visibility window: UTC \n')

for i = 1:2
    cspice_et2utc(tstart1(i),'C',1)
    cspice_et2utc(tend1(i),'C',1)
end

fprintf('   SVALBARD visibility window: UTC \n')

for i = 1:4
    cspice_et2utc(tstart2(i),'C',1)
    cspice_et2utc(tend2(i),'C',1)
end

% title('satellite propagated position from KOUROU')
figure(1)
hold on
grid on
plot(azimuth1(index_vis1),elevation1(index_vis1),'r*')
ylim([minelevation1 90])
% Plot settings
set(gca,'FontSize',12) 
legend('satellite position','Interpreter','latex','FontSize',14)
xlabel('$azimuth$ [deg]','Interpreter','latex','FontSize', 20)
ylabel('$elevation$ [deg]','Interpreter','latex','FontSize', 20)

% title('satellite propagated position from SVALBARD') 
figure(2)
hold on
grid on
plot(azimuth2(index_vis2),elevation2(index_vis2),'r*')
ylim([minelevation2 90])
% Plot settings
set(gca,'FontSize',12)
legend('satellite position','Interpreter','latex','FontSize',14)
xlabel('$azimuth$ [deg]','Interpreter','latex','FontSize', 20)
ylabel('$elevation$ [deg]','Interpreter','latex','FontSize', 20)

%% REQUEST 2

addpath('sgp4');

whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

% spacecraftName = 'Mango';
spacecraftID = 36599;

% Constant for arcseconds to radians conversions
arcsec2rad = pi / (180*3600);

% Extract the TLE from a 3LE file and initialize the inputs for SGP4
% Mango
mango = read_3LE(spacecraftID, 'tle\36599.3le', whichconst);

% Get TLE epoch
[year,mon,day,hr,min,sec] = invjday(mango.jdsatepoch, mango.jdsatepochf);
mango_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
mango_epoch_et = cspice_str2et(mango_epoch_str);

% fprintf('TLE reference epoch: UTC %s \n \n', mango_epoch_str);

ddpsi = -0.073296*arcsec2rad; %  [rad]
ddeps = -0.009373*arcsec2rad; %  [rad]

% Loop over epochs
reci_mango = zeros(3,length(tvec));
veci_mango = zeros(3,length(tvec));

for i = 1:length(tvec)

    % SGP4 propagation
    tsince = (tvec(i) - mango_epoch_et)/60.0; % minutes from TLE epoch
    [~,rteme_mango,vteme_mango] = sgp4(mango, tsince);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt = cspice_unitim(tvec(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [reci_mango(:,i), veci_mango(:,i)] = ...
        teme2eci(rteme_mango, vteme_mango, [0.0;0.0;0.0], ttt, ddpsi, ddeps);

end

% Antenna pointings
[azimuth1,elevation1,rho1,~] = antenna_pointing(station1Name, tvec,...
    [reci_mango',veci_mango']);
[azimuth2,elevation2,rho2,~] = antenna_pointing(station2Name, tvec,...
    [reci_mango',veci_mango']);

% adding noise
sigma_rho = 0.01;
sigma_angles1 = 100e-3;
sigma_angles2 = 125e-3;

sigma1 = [sigma_rho, sigma_angles1, sigma_angles1];
sigma2 = [sigma_rho, sigma_angles2, sigma_angles2];

meas1 = [rho1, azimuth1, elevation1];
meas2 = [rho2, azimuth2, elevation2];

noise_cov1 = diag(sigma1.^2);
noise_cov2 = diag(sigma2.^2);

noised_meas1 = mvnrnd(meas1,noise_cov1);
noised_meas2 = mvnrnd(meas2,noise_cov2);

rho_noise1 = noised_meas1(:,1);
az_noise1 = noised_meas1(:,2);
el_noise1 = noised_meas1(:,3);
rho_noise2 = noised_meas2(:,1);
az_noise2 = noised_meas2(:,2);
el_noise2 = noised_meas2(:,3);

% title('satellite measured position from KOUROU')
figure(3)
hold on
grid on
plot(az_noise1(index_vis1),el_noise1(index_vis1),'r*')
ylim([minelevation1 90])
% Plot settings
set(gca,'FontSize',12) 
legend('satellite position','Interpreter','latex','FontSize',14)
xlabel('$azimuth$ [deg]','Interpreter','latex','FontSize', 20)
ylabel('$elevation$ [deg]','Interpreter','latex','FontSize', 20)

% title('satellite measured position from SVALBARD') 
figure(4)
hold on
grid on
plot(az_noise2(index_vis2),el_noise2(index_vis2),'r*')
ylim([minelevation2 90])
% Plot settings
set(gca,'FontSize',12)
legend('satellite position','Interpreter','latex','FontSize',14)
xlabel('$azimuth$ [deg]','Interpreter','latex','FontSize', 20)
ylabel('$elevation$ [deg]','Interpreter','latex','FontSize', 20)

%% REQUEST 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% options
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt',...
    'Display', 'iter', 'TolFun', 1e-15, 'MaxIter',3000);  

%rng shuffle % For reproducibility
rng default % For reproducibility
ms = MultiStart;

% Cost function
fun1 = @(x) cost1(x, station1Name, tvec, index_vis1, noised_meas1, sigma1);

x0_guess = [reci_mango(:,1); veci_mango(:,1)];

problem = createOptimProblem('lsqnonlin','x0',x0_guess,...
    'objective',fun1,'lb',[],'ub',[],'options',options);

xA = run(ms,problem,15);

[xA,resnormA,residualA,exitflagA,~,~,jacA] = lsqnonlin(fun1, xA,...
    [], [], options);

% Resulting covariance
JacA= full(jacA);
P_A = resnormA / (length(residualA)-length(x0_guess)) .* inv(JacA'*JacA);
sigma_rA = sqrt(trace(P_A(1:3,1:3)));
sigma_vA = sqrt(trace(P_A(4:6,4:6)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% options
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt',...
    'Display', 'iter', 'TolFun',1e-15, 'MaxIter', 3000);  

% Cost function
fun2 = @(x) cost2(x,station1Name,station2Name,tvec,index_vis1,index_vis2,...
    noised_meas1,noised_meas2,sigma1,sigma2);

x0_guess = [reci_mango(:,1); veci_mango(:,1)];

problem = createOptimProblem('lsqnonlin','x0',x0_guess,...
    'objective',fun2,'lb',[],'ub',[],'options',options);

xB = run(ms,problem,15);

[xB,resnormB,residualB,exitflagB,~,~,jacB] = lsqnonlin(fun2, xB,...
    [], [], options);

%Resulting covariance
JacB = full(jacB);
P_B = resnormB / (length(residualB)-length(x0_guess)) .* inv(JacB'*JacB);
sigma_rB=sqrt(trace(P_B(1:3,1:3)));
sigma_vB=sqrt(trace(P_B(4:6,4:6)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% options
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt',...
    'Display', 'iter', 'TolFun', 1e-15, 'MaxIter', 3000);  

% Cost function
fun3 = @(x) cost3(x,station1Name,station2Name,tvec,index_vis1,index_vis2,...
    noised_meas1,noised_meas2,sigma1,sigma2);

x0_guess = [reci_mango(:,1); veci_mango(:,1)];

problem = createOptimProblem('lsqnonlin','x0',x0_guess,...
    'objective',fun3,'lb',[],'ub',[],'options',options);

xC = run(ms,problem,15);

[xC,resnormC,residualC,exitflagC,~,~,jacC] = lsqnonlin(fun3, xC,...
    [], [], options);

%Resulting covariance
JacC = full(jacC);
P_C = resnormC / (length(residualC)-length(x0_guess)) .* inv(JacC'*JacC);
sigma_rC=sqrt(trace(P_C(1:3,1:3)));
sigma_vC=sqrt(trace(P_C(4:6,4:6)));

% Clear kernels
cspice_kclear();

%% REQUEST 5

% Clear memory workspace and set path
clearvars; close all; clc 

% Load SPICE kernels: 
cspice_furnsh('assignment02.tm');
addpath('sgp4');
format long g

whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

t0 = '2010-08-12-05:30:00.000 UTC'; 
tf = '2010-08-12-11:00:00.000 UTC'; 
et0 = cspice_str2et(t0);  % Ephemeris time [s]
etf = cspice_str2et(tf);  % Ephemeris time [s]
etsep = cspice_str2et('2010-08-12-05:27:39.114 UTC');

mu = cspice_bodvrd('Earth','GM',1);  % Earth gravitational parameter [km^3/s^2]
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

r0_sep = [4621.69343340281; 5399.26386352847; -3.09039248714313];
v0_sep = [0.813960847513811; -0.719449862738607; 7.42706066911294];
x0_sep = [r0_sep;v0_sep];

[~, xx] = ode113(@(t,x) kepler(t,x,mu), [etsep et0], x0_sep, options);
x0 = xx(end,:);

% Define stations name 
station1Name = 'KOUROU';
station2Name = 'SVALBARD';

spacecraftName = 'Tango';
spacecraftID = 36827;

% Equally spaced vector, with one minute time-step
tvec = et0:60:etf;

[time, tangostate] = ode113(@(t,x) kepler(t,x,mu), tvec, x0, options);

[~,elevation1,~,~] = antenna_pointing(station1Name, tvec, tangostate);
[~,elevation2,~,~] = antenna_pointing(station2Name, tvec, tangostate);

minelevation1 = 10; % deg
minelevation2 = 5; % deg

index_vis1 = elevation1 >= minelevation1; % KOUROU visibility window
index_vis2 = elevation2 >= minelevation2; % SVALBARD visibility window

% Constant for arcseconds to radians conversions
arcsec2rad = pi / (180*3600);

% Extract the TLE from a 3LE file and initialize the inputs for SGP4
% Tango
tango = read_3LE(spacecraftID, 'tle\36827.3le', whichconst);

% Get TLE epoch
[year,mon,day,hr,min,sec] = invjday(tango.jdsatepoch, tango.jdsatepochf);
tango_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
tango_epoch_et = cspice_str2et(tango_epoch_str);

ddpsi = -0.073296*arcsec2rad; %  [rad]
ddeps = -0.009373*arcsec2rad; %  [rad]

% Loop over epochs
reci_tango = zeros(3,length(tvec));
veci_tango = zeros(3,length(tvec));

for i = 1:length(tvec)

    % SGP4 propagation
    tsince = (tvec(i) - tango_epoch_et)/60.0; % minutes from TLE epoch
    [~,rteme_tango,vteme_tango] = sgp4(tango, tsince);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt = cspice_unitim(tvec(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [reci_tango(:,i), veci_tango(:,i)] = ...
        teme2eci(rteme_tango, vteme_tango, [0.0;0.0;0.0], ttt, ddpsi, ddeps);

end

% Antenna pointings

[azimuth1,elevation1,rho1,~] = antenna_pointing(station1Name, tvec,...
    [reci_tango',veci_tango']);
[azimuth2,elevation2,rho2,~] = antenna_pointing(station2Name, tvec,...
    [reci_tango',veci_tango']);

% adding noise
sigma_rho = 0.01;
sigma_angles1 = 100e-3;
sigma_angles2 = 125e-3;

sigma1 = [sigma_rho, sigma_angles1, sigma_angles1];
sigma2 = [sigma_rho, sigma_angles2, sigma_angles2];

meas1 = [rho1, azimuth1, elevation1];
meas2 = [rho2, azimuth2, elevation2];

noise_cov1 = diag(sigma1.^2);
noise_cov2 = diag(sigma2.^2);

noised_meas1 = mvnrnd(meas1,noise_cov1);
noised_meas2 = mvnrnd(meas2,noise_cov2);

% options
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt',...
    'Display', 'iter', 'TolFun', 1e-15, 'MaxIter', 3000);  

rng default % For reproducibility
ms = MultiStart;

% Cost function
fun = @(x) cost3(x,station1Name,station2Name,tvec,index_vis1,index_vis2,...
    noised_meas1,noised_meas2,sigma1,sigma2);

x0_guess = [reci_tango(:,1); veci_tango(:,1)];

problem = createOptimProblem('lsqnonlin','x0',x0_guess,...
    'objective',fun,'lb',[],'ub',[],'options',options);

x = run(ms,problem,15);

[x,resnorm,residual,exitflag,~,~,jac] = lsqnonlin(fun, x,...
    [], [], options);

%Resulting covariance
Jac = full(jac);
P = resnorm / (length(residual)-length(x0_guess)) .* inv(Jac'*Jac);
sigma_r = sqrt(trace(P(1:3,1:3)));
sigma_v = sqrt(trace(P(4:6,4:6)));

% Clear kernels
cspice_kclear();


%% FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dxdt] = kepler(~, state, mu)
%KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian) propagator
%   Evaluates the right-hand-side of a newtonian 2-body propagator.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 24/10/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   t   : [ 1, 1] epoch (unused)
%   x   : [6, 1] cartesian state vector wrt Solar-System-Barycentre and
%                 State Transition Matrix elements
%   GM  : [ 1, 1] gravitational constant of the body
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
%

% Initialize right-hand-side
dxdt = zeros(6,1);

x = state(1);
y = state(2);
z = state(3);
rr = [x; y; z]; % column

% Compute distance
dist = sqrt(x^2 + y^2 + z^2);

% Position detivative is object's velocity
dxdt(1:3) = state(4:6);   
% Compute the gravitational acceleration using Newton's law
dxdt(4:6) = - mu * rr /(dist^3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dxdt] = keplerJ2(t, state, mu)
% The function evaluates the right-hand-side of a 2-body (keplerian)
% propagator considering also the J2 perturbation effect
%   
% Author: Alessandro Cambielli
%  
% Inputs:
%   t : time [1,1]
%   state : state of the body (position and velocity) [6,1]
%   mu : gravitational parameter [1,1]    
%
% Outputs:
%   dxdt : [6,1] 
%

% Initialize right-hand-side
dxdt = zeros(6,1);

x = state(1);
y = state(2);
z = state(3);
rr = [x; y; z]; % ECI position

J2 = 0.0010826269;
Re = cspice_bodvrd('Earth','RADII',3);

ECI2ECEF = cspice_pxform('J2000', 'ITRF93', t);
ECEF2ECI = cspice_pxform('ITRF93', 'J2000', t);

rr_ECEF = ECI2ECEF*rr; % ECEF position

% Compute distance ECEF
dist_ECEF = sqrt(rr_ECEF(1)^2 + rr_ECEF(2)^2 + rr_ECEF(3)^2);

% Compute distance ECI
dist_ECI = sqrt(x^2 + y^2 + z^2);

a_J2_ECEF = 3/2*mu*J2/dist_ECEF^3*rr_ECEF*(Re(1)/dist_ECEF)^2.*...
    (5*(rr_ECEF(3)/dist_ECEF)^2-[1;1;3]); % J2 ECEF
a_J2_ECI = ECEF2ECI*a_J2_ECEF;

% Position detivative is object's velocity
dxdt(1:3) = state(4:6);

% Compute the gravitational acceleration using Newton's law
dxdt(4:6) = - mu * rr /(dist_ECI^3) + a_J2_ECI;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [azimuth,elevation,range,rangerate] = antenna_pointing(stationName,...
    tvec, mangostate)
% The function evaluates the azimuth, elevation, range and range-rate of a
% body as perceived from a station in a certain time window
%   
% Author: Alessandro Cambielli
%  
% Inputs:
%   statioName : name of the considered station [string]
%   tvec : time interval of the observation [1,n]
%   mangostate : state of the satellite in ECI [n,6]    
%
% Outputs:
%   azimuth : azimuth of the satellite [n,1] 
%   elevation : elevation of the satellite [n,1]
%   range : range of the satellite [n,1]
%   range-rate : range-rate of the satellite [n,1]
%

topoFrame = [stationName, '_TOPO'];

range = zeros(length(tvec),1);
azimuth = zeros(length(tvec),1);
elevation = zeros(length(tvec),1);
rangerate = zeros(length(tvec),1);

for i = 1:length(tvec)

    % ECI2TOPO matrix using pinpoint
    ECI2TOPO = cspice_sxform('J2000', topoFrame, tvec(i));

    % Station position in ECI
    antenna_ECI = cspice_spkezr(stationName, tvec(i), 'J2000',...
    'NONE', 'EARTH');
    
    % relative position in ECI and TOPO
    relative_ECI = mangostate(i,:)' - antenna_ECI;
    relative_TOPO = ECI2TOPO*relative_ECI;

    % retrieve parameters 
    parameters = cspice_xfmsta(relative_TOPO,...
    'RECTANGULAR','LATITUDINAL','EARTH');
    range(i) = parameters(1);
    azimuth(i) = rad2deg(parameters(2));
    elevation(i) = rad2deg(parameters(3));
    rangerate(i) = parameters(4);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function residual = cost1(x, stationName, tvec, window, meas_real, sigma)
% Cost function considering a single station and keplerian motion
%   
% Author: Alessandro Cambielli
%  
% Inputs:
%   x : initial state of the satellite [6,1]
%   statioName : name of the considered station [string]
%   tvec : time interval of the observation [1,n]
%   window : visibility window for the station [n,1]
%   meas_real : noised measurements from the station [n,3]
%   sigma : measurements noises [1,3]    
%
% Outputs:
%   residual : residual for the solver [1,1] 
%

mu = cspice_bodvrd('Earth','GM',1); 
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

[~,mangostate] = ode113(@(t,xx) kepler(t,xx,mu),tvec,x,options);

[azimuth,elevation,rho,~] = antenna_pointing(stationName,tvec,mangostate);
meas_pred(:,1) = rho(window);
meas_pred(:,2) = azimuth(window);
meas_pred(:,3) = elevation(window);

rho_real = meas_real(:,1);
rho_real = rho_real(window);
az_real = meas_real(:,2);
az_real = az_real(window);
el_real = meas_real(:,3);
el_real = el_real(window);
meas_real = [rho_real, az_real, el_real];

% Compute the residual
W_m = diag(1./sigma);
weighted_res = W_m*(meas_pred - meas_real)';
residual = weighted_res(:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function residual = cost2(x,stationName1,stationName2,tvec,window1,window2,...
    meas_real1,meas_real2,sigma1,sigma2)
% Cost function considering two separate stations and keplerian motion
%   
% Author: Alessandro Cambielli
%  
% Inputs:
%   x : initial state of the satellite [6,1]
%   statioName1 : name of the first considered station [string]
%   statioName2 : name of the second considered station [string]
%   tvec : time interval of the observation [1,n]
%   window1 : visibility window for the first station [n,1]
%   window2 : visibility window for the second station [n,1]
%   meas_real1 : noised measurements from the first station [n,3]
%   meas_real2 : noised measurements from the second station [n,3]
%   sigma1 : first station measurements noises [1,3]
%   sigma2 : second station measurements noises [1,3]
%
% Outputs:
%   residual : residual for the solver [1,1] 
%

mu = cspice_bodvrd('Earth','GM',1); 
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

[~,mangostate] = ode113(@(t,xx) kepler(t,xx,mu),tvec,x,options);

[azimuth1,elevation1,rho1,~] = antenna_pointing(stationName1,tvec,mangostate);
meas_pred1(:,1) = rho1(window1);
meas_pred1(:,2) = azimuth1(window1);
meas_pred1(:,3) = elevation1(window1);

[azimuth2,elevation2,rho2,~] = antenna_pointing(stationName2,tvec,mangostate);
meas_pred2(:,1) = rho2(window2);
meas_pred2(:,2) = azimuth2(window2);
meas_pred2(:,3) = elevation2(window2);

rho_real1 = meas_real1(:,1);
rho_real1 = rho_real1(window1);
az_real1 = meas_real1(:,2);
az_real1 = az_real1(window1);
el_real1 = meas_real1(:,3);
el_real1 = el_real1(window1);
meas_real1 = [rho_real1, az_real1, el_real1];

rho_real2 = meas_real2(:,1);
rho_real2 = rho_real2(window2);
az_real2 = meas_real2(:,2);
az_real2 = az_real2(window2);
el_real2 = meas_real2(:,3);
el_real2 = el_real2(window2);
meas_real2 = [rho_real2, az_real2, el_real2];

% Compute the residual
W_m1 = diag(1./sigma1);
weighted_res1 = W_m1*(meas_pred1 - meas_real1)';

W_m2 = diag(1./sigma2);
weighted_res2 = W_m2*(meas_pred2 - meas_real2)';

residual = [weighted_res1(:); weighted_res2(:)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function residual = cost3(x,stationName1,stationName2,tvec,window1,window2,...
    meas_real1,meas_real2,sigma1,sigma2)
% Cost function considering two separate stations and keplerian J2 
% perturbed motion   
%
% Author: Alessandro Cambielli
%  
% Inputs:
%   x : initial state of the satellite [6,1]
%   statioName1 : name of the first considered station [string]
%   statioName2 : name of the second considered station [string]
%   tvec : time interval of the observation [1,n]
%   window1 : visibility window for the first station [n,1]
%   window2 : visibility window for the second station [n,1]
%   meas_real1 : noised measurements from the first station [n,3]
%   meas_real2 : noised measurements from the second station [n,3]
%   sigma1 : first station measurements noises [1,3]
%   sigma2 : second station measurements noises [1,3]
%
% Outputs:
%   residual : residual for the solver [1,1] 
%

mu = cspice_bodvrd('Earth','GM',1); 
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

[~,mangostate] = ode113(@(t,xx) keplerJ2(t,xx,mu),tvec,x,options);

[azimuth1,elevation1,rho1,~] = antenna_pointing(stationName1,tvec,mangostate);
meas_pred1(:,1) = rho1(window1);
meas_pred1(:,2) = azimuth1(window1);
meas_pred1(:,3) = elevation1(window1);

[azimuth2,elevation2,rho2,~] = antenna_pointing(stationName2,tvec,mangostate);
meas_pred2(:,1) = rho2(window2);
meas_pred2(:,2) = azimuth2(window2);
meas_pred2(:,3) = elevation2(window2);

rho_real1 = meas_real1(:,1);
rho_real1 = rho_real1(window1);
az_real1 = meas_real1(:,2);
az_real1 = az_real1(window1);
el_real1 = meas_real1(:,3);
el_real1 = el_real1(window1);
meas_real1 = [rho_real1, az_real1, el_real1];

rho_real2 = meas_real2(:,1);
rho_real2 = rho_real2(window2);
az_real2 = meas_real2(:,2);
az_real2 = az_real2(window2);
el_real2 = meas_real2(:,3);
el_real2 = el_real2(window2);
meas_real2 = [rho_real2, az_real2, el_real2];

% Compute the residual
W_m1 = diag(1./sigma1);
weighted_res1 = W_m1*(meas_pred1 - meas_real1)';

W_m2 = diag(1./sigma2);
weighted_res2 = W_m2*(meas_pred2 - meas_real2)';

residual = [weighted_res1(:); weighted_res2(:)];

end


