%% Spacecraft Guidance and Navigation - Periodic Orbit (2023/2024)
% Assignment:     2
% Exercise:       3 - Sequential filters
% Author:         Alessandro Cambielli

% ID number:      10619322
%                 102615

%% REQUEST 1

% Clear memory workspace and set path
clearvars; close all; clc 

% Load SPICE kernels: 
cspice_furnsh('assignment02.tm');
addpath('sgp4');
format long g

whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

t0 = '2010-08-12-05:30:00.000 UTC'; 
tf = '2010-08-12-06:30:00.000 UTC'; 
et0 = cspice_str2et(t0);  % Ephemeris time [s]
etf = cspice_str2et(tf);  % Ephemeris time [s]
etsep = cspice_str2et('2010-08-12-05:27:39.114 UTC'); % Ephemeris time [s]

mu = cspice_bodvrd('Earth','GM',1);  % Earth gravitational parameter [km^3/s^2]
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

r0_sep = [4622.232026629; 5399.3369588058; -0.0212138165769957];
v0_sep = [0.812221125483763; -0.721512914578826; 7.42665302729053];
x0_sep = [r0_sep;v0_sep];

P01 = [5.6e-7 3.5e-7 -7.1e-8 0 0 0;...
      3.5e-7 9.7e-7 7.6e-8 0 0 0;...
      -7.1e-8 7.6e-8 8.1e-8 0 0 0;...
      0 0 0 2.8e-11 0 0;...
       0 0 0 0 2.7e-11 0;...
       0 0 0 0 0 9.6e-12];

[~, xx] = ode113(@(t,x) keplerJ2(t,x,mu), [etsep et0], x0_sep, options);
x0 = xx(end,:);

% Define stations name 
stationName = 'SVALBARD';

%spacecraftName = 'Mango';
spacecraftID1 = 36599;
%spacecraftName = 'Tango';
spacecraftID2 = 36827;

% Equally spaced vector, with 5 seconds time-step
tvec = et0:5:etf;

[~, mangostate] = ode113(@(t,x) keplerJ2(t,x,mu), tvec, x0, options);

[~,elevation1,~,~] = antenna_pointing(stationName, tvec, mangostate);

minelevation = deg2rad(5); % deg

index_vis = elevation1 >= minelevation; % SVALBARD visibility window

tstart = tvec(168);
tend = tvec(296);

fprintf('   SVALBARD visibility window: UTC \n')

cspice_et2utc(tstart,'C',1)
cspice_et2utc(tend,'C',1)

% Constant for arcseconds to radians conversions
arcsec2rad = pi / (180*3600);

% Extract the TLE from a 3LE file and initialize the inputs for SGP4
% Mango
mango = read_3LE(spacecraftID1, 'tle\36599.3le', whichconst);

% Get TLE epoch
[year,mon,day,hr,min,sec] = invjday(mango.jdsatepoch, mango.jdsatepochf);
mango_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
mango_epoch_et = cspice_str2et(mango_epoch_str);

ddpsi = -0.073296*arcsec2rad; %  [rad]
ddeps = -0.009373*arcsec2rad; %  [rad]

tvec1 = tvec(index_vis);

% Loop over epochs
reci_mango = zeros(3,length(tvec1));
veci_mango = zeros(3,length(tvec1));

for i = 1:length(tvec1)

    % SGP4 propagation
    tsince = (tvec1(i) - mango_epoch_et)/60.0; % minutes from TLE epoch
    [~,rteme_mango,vteme_mango] = sgp4(mango, tsince);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt = cspice_unitim(tvec1(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [reci_mango(:,i), veci_mango(:,i)] = ...
        teme2eci(rteme_mango, vteme_mango, [0.0;0.0;0.0], ttt, ddpsi, ddeps);

end

% Antenna pointings

[azimuth_TLE,elevation_TLE,rho_TLE,~] = antenna_pointing(stationName, tvec1,...
    [reci_mango',veci_mango']);

% adding noise
sigma_rho1 = 0.01; % km
sigma_angles1 = deg2rad(125e-3); % rad

sigma1 = [sigma_angles1, sigma_angles1, sigma_rho1];
meas1 = [azimuth_TLE, elevation_TLE, rho_TLE];

noise_cov1 = diag(sigma1.^2);
noised_meas1 = mvnrnd(meas1,noise_cov1);

% Unscented Kalman Filter

R1 = diag([sigma_angles1^2 sigma_angles1^2 sigma_rho1^2]);

mean1 = zeros(6,length(tvec1));
P1 = zeros(6,6,length(tvec1));

[mean1(:,1),P1(:,:,1)] = UT(x0_sep,P01*1e4,etsep,tvec1(1));

for i = 2:length(tvec1)

    [mean1(:,i),P1(:,:,i)] = UKF(mean1(:,i-1),P1(:,:,i-1),...
        tvec1(i-1),tvec1(i),stationName,R1,noised_meas1(i,:)');

end

% compute errors
pos_error1 = zeros(1,length(tvec1));
vel_error1 = zeros(1,length(tvec1));
pos_3sigma1 = zeros(1,length(tvec1));
vel_3sigma1 = zeros(1,length(tvec1));

for i = 1:length(tvec1)
    pos_error1(i) = norm(reci_mango(:,i) - mean1(1:3,i));
    vel_error1(i) = norm(veci_mango(:,i) - mean1(4:6,i));
    pos_3sigma1(i) = 3*sqrt(trace(P1(1:3,1:3,i)));
    vel_3sigma1(i) = 3*sqrt(trace(P1(4:6,4:6,i)));
end

% plots

epoch = linspace(datetime(2010,08,12,05,43,55),datetime(2010,08,12,05,54,35),...
    length(tvec1));

% title('Mango position error and 3 sigma')
figure(1)
plot(epoch,pos_error1,'r-')
hold on
grid on
plot(epoch,pos_3sigma1,'b-')
xlim([epoch(1) epoch(end)])
% Plot settings
set(gca,'FontSize',12) 
legend('position error','$3\sigma$','Interpreter','latex','FontSize',14)
xlabel('$epoch$ [date]','Interpreter','latex','FontSize', 20)
ylabel('$error$ [km]','Interpreter','latex','FontSize', 20)

% title('Mango velocity error and 3 sigma')
figure(2)
plot(epoch,vel_error1,'r-')
hold on
grid on
plot(epoch,vel_3sigma1,'b-')
xlim([epoch(1) epoch(end)])
% Plot settings
set(gca,'FontSize',12) 
legend('velocity error','$3\sigma$','Interpreter','latex','FontSize',14)
xlabel('$epoch$ [date]','Interpreter','latex','FontSize', 20)
ylabel('$error$ [km/s]','Interpreter','latex','FontSize', 20)

%% REQUEST 2

% Extract the TLE from a 3LE file and initialize the inputs for SGP4
% Mango
mango = read_3LE(spacecraftID1, 'tle\36599.3le', whichconst);
% Tango
tango = read_3LE(spacecraftID2, 'tle\36827.3le', whichconst);

% Get TLE epoch for Mango
[year,mon,day,hr,min,sec] = invjday(mango.jdsatepoch, mango.jdsatepochf);
mango_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
mango_epoch_et = cspice_str2et(mango_epoch_str);

% Get TLE epoch for Tango
[year,mon,day,hr,min,sec] = invjday(tango.jdsatepoch, tango.jdsatepochf);
tango_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
tango_epoch_et = cspice_str2et(tango_epoch_str);

% Equally spaced vector, with 1 second time-step
tvecmango = mango_epoch_et:1:et0;
tvectango = tango_epoch_et:1:et0;

% Loop over epochs
reci_mango = zeros(3,length(tvecmango));
veci_mango = zeros(3,length(tvecmango));
reci_tango = zeros(3,length(tvectango));
veci_tango = zeros(3,length(tvectango));

for i = 1:length(tvecmango)

    % SGP4 propagation
    tsince_mango = (tvecmango(i) - mango_epoch_et)/60.0; % minutes from TLE epoch
    [~,rteme_mango,vteme_mango] = sgp4(mango, tsince_mango);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    tttmango = cspice_unitim(tvecmango(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [reci_mango(:,i), veci_mango(:,i)] = ...
        teme2eci(rteme_mango, vteme_mango, [0.0;0.0;0.0], tttmango, ddpsi, ddeps);

end

for i = 1:length(tvectango)

    % SGP4 propagation
    tsince_tango = (tvectango(i) - tango_epoch_et)/60.0;
    [~,rteme_tango,vteme_tango] = sgp4(tango, tsince_tango);

    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttttango = cspice_unitim(tvectango(i), 'ET', 'TDT')/cspice_jyear()/100;

    % TEME to ECI conversion
    [reci_tango(:,i), veci_tango(:,i)] = ...
        teme2eci(rteme_tango, vteme_tango, [0.0;0.0;0.0], ttttango, ddpsi, ddeps);

end

relative_posECI =  reci_tango(:,end) - reci_mango(:,end);
relative_velECI =  veci_tango(:,end) - veci_mango(:,end);

% LVLH frame
r_ECI = reci_mango(:,end);
r_hat_ECI = r_ECI/norm(r_ECI);
v_ECI = veci_mango(:,end);
v_hat_ECI = v_ECI/norm(v_ECI);
h_ECI = cross(r_ECI,v_ECI);
h_hat_ECI = h_ECI/norm(h_ECI);

i = r_hat_ECI;
k = h_hat_ECI;
j = cross(k,i);

idot = 1/norm(r_ECI)*(v_ECI-dot(r_hat_ECI,v_ECI)*r_hat_ECI);
kdot = zeros(3,1);
jdot = cross(-h_hat_ECI,idot);

ECI2LVLH = [i, j, k]';
ECI2LVLHdot = [idot, jdot, kdot]';

posLVLH = ECI2LVLH*relative_posECI;
velLVLH = ECI2LVLHdot*relative_posECI + ECI2LVLH*relative_velECI;

x0_rel = [posLVLH;velLVLH];

R = norm(r_ECI);
n = sqrt(mu/(R^3));

% Equally spaced vector, with 5 seconds time-step
tvec = et0:5:etf;

[~,relstate] = ode113(@(t,x) CW(t,x,n),tvec,x0_rel,options);

[azimuth2,elevation2,range2] = FFRF(tvec,relstate);

% adding noise
sigma_angles2 = deg2rad(1); % rad
sigma_rho2 = 1e-5; % km

sigma2 = [sigma_angles2, sigma_angles2, sigma_rho2];
meas2 = [azimuth2, elevation2, range2];

noise_cov2 = diag(sigma2.^2);
noised_meas2 = mvnrnd(meas2,noise_cov2);

P02 = diag([0.01, 0.01, 0.1, 0.0001, 0.0001, 0.001]);

tvec2 = tvec1(end)+5:5:(tvec1(end)+5+20*60);
window = find(index_vis);
index = window(end);

R2 = diag([sigma_angles2^2 sigma_angles2^2 sigma_rho2^2]);

mean2 = zeros(6,length(tvec2));
P2 = zeros(6,6,length(tvec2));

mean2(:,1) = relstate(index+1,:);
P2(:,:,1) = P02;

for i = 2:length(tvec2)

    [mean2(:,i),P2(:,:,i)] = UKF2(n,mean2(:,i-1),P2(:,:,i-1),...
        tvec2(i-1),tvec2(i),R2,noised_meas2(i+index-1,:)');

end

% compute errors
pos_error2 = zeros(1,length(tvec2));
vel_error2 = zeros(1,length(tvec2));
pos_3sigma2 = zeros(1,length(tvec2));
vel_3sigma2 = zeros(1,length(tvec2));

for i = 1:length(tvec2)
    pos_error2(i) = norm(relstate(i+index-1,1:3)' - mean2(1:3,i));
    vel_error2(i) = norm(relstate(i+index-1,4:6)' - mean2(4:6,i));
    pos_3sigma2(i) = 3*sqrt(trace(P2(1:3,1:3,i)));
    vel_3sigma2(i) = 3*sqrt(trace(P2(4:6,4:6,i)));
end

% plots

epoch = linspace(datetime(2010,08,12,05,54,35),datetime(2010,08,12,06,14,35),...
    length(tvec2));

% title('relative position error and 3 sigma')
figure(3)
plot(epoch,pos_error2,'r-')
hold on
grid on
plot(epoch,pos_3sigma2,'b-')
xlim([epoch(1) epoch(end)])
% Plot settings
set(gca,'FontSize',12) 
legend('position error','$3\sigma$','Interpreter','latex','FontSize',14)
xlabel('$epoch$ [date]','Interpreter','latex','FontSize', 20)
ylabel('$error$ [km]','Interpreter','latex','FontSize', 20)

% title('relative velocity error and 3 sigma')
figure(4)
plot(epoch,vel_error2,'r-')
hold on
grid on
plot(epoch,vel_3sigma2,'b-')
xlim([epoch(1) epoch(end)])
ylim([-0.001 0.12])
% Plot settings
set(gca,'FontSize',12) 
legend('velocity error','$3\sigma$','Interpreter','latex','FontSize',14)
xlabel('$epoch$ [date]','Interpreter','latex','FontSize', 20)
ylabel('$error$ [km/s]','Interpreter','latex','FontSize', 20)

%% REQUEST 3

meanTango1 = zeros(length(tvec2),6);
PTango1 = zeros(6,6,length(tvec2));

meanTango1(1,:) = mean1(:,end)';
PTango1(:,:,1) = P1(:,:,end);

for i = 2:length(tvec2)

    [meanTango1(i,:),PTango1(:,:,i)] = UT(meanTango1(i-1,:)',...
        PTango1(:,:,i-1),tvec2(i-1),tvec2(i));
   
end

LVLH2ECI = ECI2LVLH';
LVLH2ECIdot = ECI2LVLHdot';
psi = [LVLH2ECI, zeros(3); LVLH2ECIdot, LVLH2ECI];

PTango2 = zeros(6,6,length(tvec2));
Ptango = zeros(6,6,length(tvec2));

for i = 1:length(tvec2)
    PTango2(:,:,i) = psi*P2(:,:,i)*psi';
end

for i = 1:length(tvec2)
    Ptango(:,:,i) = PTango1(:,:,i) + PTango2(:,:,i);
end

pos_3sigma3 = zeros(1,length(tvec2));
vel_3sigma3 = zeros(1,length(tvec2));

for i = 1:length(tvec2)
    pos_3sigma3(i) = 3*sqrt(trace(Ptango(1:3,1:3,i)));
    vel_3sigma3(i) = 3*sqrt(trace(Ptango(4:6,4:6,i)));
end

% title('Tango position 3 sigma')
figure(5)
plot(epoch(2:end),pos_3sigma3(2:end),'b-')
grid on
xlim([epoch(1) epoch(end)])
% Plot settings
set(gca,'FontSize',12) 
legend('$position$ $3\sigma$','Interpreter','latex','FontSize',14)
xlabel('$epoch$ [date]','Interpreter','latex','FontSize', 20)
ylabel('$3\sigma$ [km]','Interpreter','latex','FontSize', 20)

% title('Tango velocity 3 sigma')
figure(6)
plot(epoch(2:end),vel_3sigma3(2:end),'b-')
grid on
xlim([epoch(1) epoch(end)])
% Plot settings
set(gca,'FontSize',12) 
legend('$velocity$ $3\sigma$','Interpreter','latex','FontSize',14)
xlabel('$epoch$ [date]','Interpreter','latex','FontSize', 20)
ylabel('$3\sigma$ [km/s]','Interpreter','latex','FontSize', 20)


%% FUNCTIONS

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
%   dxdt : derivative of the state vector [6,1] 
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

function [mean,P] = UT(mean0,P0,t0,tf)
% The function propagates initial mean state and covariance exploiting the
% unscented transform.
%
% Author: Alessandro Cambielli
%
% Inputs:
%   mean0 : initial mean state [6,1] 
%   P0 : initial covariance matrix [6,6] 
%   t0 : initial propagation time [1,1] 
%   tf : final propagation time [1,1]
%
% Outputs:
%   mean : final mean state [6,1] 
%   P : final covariance matrix [6,6]
%

mu = cspice_bodvrd('Earth','GM',1);
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

n = 6;
alpha = 1;
k = 0;
lambda = alpha^2*(n+k) - n;
matrix = sqrtm((n+lambda)*P0); % 6,6

chi = zeros(6,2*n+1); % 6,13

chi(:,1) = mean0;
for i = 1:n
    chi(:,i+1) = mean0 + matrix(:,i);
    chi(:,i+1+n) = mean0 - matrix(:,i);
end

beta = 2;
w_m = zeros(1,2*n+1); % 1,13
w_c = zeros(1,2*n+1); % 1,13

for i = 1:2*n+1
    if i == 1
        w_m(i) = lambda/(n+lambda);
        w_c(i) = lambda/(n+lambda) + (1-alpha^2+beta);
    else 
        w_m(i) = 1/(2*(n+lambda));
        w_c(i) = 1/(2*(n+lambda));
    end
end

mean = zeros(n,1);
P = zeros(n,n);
jota = zeros(n,2*n+1);

for i = 1:(2*n+1)
    [~, xx] = ode113(@(t,x) keplerJ2(t,x,mu), [t0 tf], chi(:,i), options);
    jota(:,i) = xx(end,:)';
    mean = mean + w_m(i) * jota(:,i);
end

for k = 1:(2*n+1)
    P = P + w_c(k) * (jota(:,k)-mean) * (jota(:,k)-mean)';
end

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
    azimuth(i) = parameters(2);
    elevation(i) = parameters(3);
    rangerate(i) = parameters(4);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [azimuth,elevation,range] = FFRF(tvec,relstate)
% The function evaluates the azimuth, elevation and range of a body as
% perceived from the FFRF in a certain time window
%   
% Author: Alessandro Cambielli
%  
% Inputs:
%   tvec : time interval of the observation [1,n]
%   relstate : relative state of the satellites in ECI [n,6]    
%
% Outputs:
%   azimuth : azimuth of the satellite [n,1] 
%   elevation : elevation of the satellite [n,1]
%   range : range of the satellite [n,1]
%

rel_xy = zeros(length(tvec),1);
range = zeros(length(tvec),1);
azimuth = zeros(length(tvec),1);
elevation = zeros(length(tvec),1);

for i = 1:length(tvec)
    % exploiting geometrical relationships
    rel_xy(i) = sqrt(relstate(i,1)^2 + relstate(i,2)^2);
    azimuth(i) = atan(relstate(i,2)/relstate(i,1));
    elevation(i) = atan(relstate(i,3)/rel_xy(i));
    range(i) = norm(relstate(i,1:3));
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [meanf,Pf] = UKF(mean0,P0,t0,tf,stationName,R,y_real)
% The function simulates the behaviour of an Unscented Kalman Filter.
%
% Author: Alessandro Cambielli
%
% Inputs:
%   mean0 : initial mean state [6,1]  
%   P0 : initial covariance matrix [6,6] 
%   t0 : initial time instant [1,1] 
%   tf : final time instant [1,1] 
%   stationName : name of the observing station [string] 
%   R : measurement error matrix [3,3] 
%   y_real : measurements from the station [3,1] 
%
% Outputs:
%   meanf : final mean state [6,1] 
%   Pf : final covariance matrix [6,6] 
%

mu = cspice_bodvrd('Earth','GM',1);
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

n = 6;
alpha = 1;
k = 0;
lambda = alpha^2*(n+k) - n;
matrix = sqrtm((n+lambda)*P0); % 6,6

chi = zeros(6,2*n+1); % 6,13

chi(:,1) = mean0;
for i = 1:n
    chi(:,i+1) = mean0 + matrix(:,i);
    chi(:,i+1+n) = mean0 - matrix(:,i);
end

beta = 2;
w_m = zeros(1,2*n+1); % 1,13
w_c = zeros(1,2*n+1); % 1,13

for i = 1:2*n+1
    if i == 1
        w_m(i) = lambda/(n+lambda);
        w_c(i) = lambda/(n+lambda) + (1-alpha^2+beta);
    else 
        w_m(i) = 1/(2*(n+lambda));
        w_c(i) = 1/(2*(n+lambda));
    end
end

% Perform integration
jota = zeros(n,2*n+1);

for i = 1:(2*n+1)
    [~, xx] = ode113(@(t,x) keplerJ2(t,x,mu), [t0 tf], chi(:,i), options);
    jota(:,i) = xx(end,:)';
end

gamma = zeros(3,2*n+1);

for i = 1:(2*n+1)
    [azimuth,elevation,range,~] = antenna_pointing(stationName,...
        tf,jota(:,i)');
    gamma(:,i) = [azimuth;elevation;range];
end

meanx = zeros(6,1);

for i = 1:(2*n+1)
    meanx = meanx + w_m(i)*jota(:,i);
end

meany = zeros(3,1);

for i = 1:(2*n+1)
    meany = meany + w_m(i)*gamma(:,i);
end

P = zeros(6);

for i = 1:(2*n+1)
    P = P + w_c(i)*(jota(:,i)-meanx)*(jota(:,i)-meanx)';
end

Pee = zeros(3);

for i = 1:(2*n+1)
    diff = angdiff(gamma(1:2,i),meany(1:2));
    diff2 = [diff; gamma(3,i)-meany(3)];
    Pee = Pee + w_c(i)*(diff2*diff2');
end 

Pee = Pee + R; 

Pxy = zeros(6,3);

for i = 1:(2*n+1)
    diff = angdiff(gamma(1:2,i),meany(1:2));
    diff2 = [diff; gamma(3,i)-meany(3)];
    Pxy = Pxy + w_c(i)*(jota(:,i)-meanx)*diff2';
end

deltaangles = angdiff(y_real(1:2),meany(1:2));
deltay = [deltaangles; y_real(3)-meany(3)];

K = Pxy/Pee;
meanf = meanx + K*deltay;
Pf = P - K*Pee*K';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [meanf,Pf] = UKF2(N,mean0,P0,t0,tf,R,y_real)
% The function simulates the behaviour of an Unscented Kalman Filter.
%
% Author: Alessandro Cambielli
%
% Inputs:
%   N : mean motion of Mango [1,1]
%   mean0 : initial mean state [6,1]  
%   P0 : initial covariance matrix [6,6] 
%   t0 : initial time instant [1,1] 
%   tf : final time instant [1,1]  
%   R : measurement error matrix [3,3] 
%   y_real : measurements from the station [3,1] 
%
% Outputs:
%   meanf : final mean state [6,1] 
%   Pf : final covariance matrix [6,6] 
%

options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

n = 6;
alpha = 1;
k = 0;
lambda = alpha^2*(n+k) - n;
matrix = sqrtm((n+lambda)*P0); % 6,6

chi = zeros(6,2*n+1); % 6,13

chi(:,1) = mean0;
for i = 1:n
    chi(:,i+1) = mean0 + matrix(:,i);
    chi(:,i+1+n) = mean0 - matrix(:,i);
end

beta = 2;
w_m = zeros(1,2*n+1); % 1,13
w_c = zeros(1,2*n+1); % 1,13

for i = 1:2*n+1
    if i == 1
        w_m(i) = lambda/(n+lambda);
        w_c(i) = lambda/(n+lambda) + (1-alpha^2+beta);
    else 
        w_m(i) = 1/(2*(n+lambda));
        w_c(i) = 1/(2*(n+lambda));
    end
end

% Perform integration
jota = zeros(n,2*n+1);

for i = 1:(2*n+1)
    [~, xx] = ode113(@(t,x) CW(t,x,N), [t0 tf], chi(:,i), options);
    jota(:,i) = xx(end,:)';
end

gamma = zeros(3,2*n+1);

for i = 1:(2*n+1)
    [azimuth,elevation,range] = FFRF(tf,jota(:,i)'); % PROBLEMA
    gamma(:,i) = [azimuth;elevation;range];
end

meanx = zeros(6,1);

for i = 1:(2*n+1)
    meanx = meanx + w_m(i)*jota(:,i);
end

meany = zeros(3,1);

for i = 1:(2*n+1)
    meany = meany + w_m(i)*gamma(:,i);
end

P = zeros(6);

for i = 1:(2*n+1)
    P = P + w_c(i)*(jota(:,i)-meanx)*(jota(:,i)-meanx)';
end

Pee = zeros(3);

for i = 1:(2*n+1)
    diff = angdiff(gamma(1:2,i),meany(1:2));
    diff2 = [diff; gamma(3,i)-meany(3)];
    Pee = Pee + w_c(i)*(diff2*diff2');
end 

Pee = Pee + R; 

Pxy = zeros(6,3);

for i = 1:(2*n+1)
    diff = angdiff(gamma(1:2,i),meany(1:2));
    diff2 = [diff; gamma(3,i)-meany(3)];
    Pxy = Pxy + w_c(i)*(jota(:,i)-meanx)*diff2';
end

deltaangles = angdiff(y_real(1:2),meany(1:2));
deltay = [deltaangles; y_real(3)-meany(3)];

K = Pxy/Pee;
meanf = meanx + K*deltay;
Pf = P - K*Pee*K';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dxdt] = CW(~,state,n)
% The function evaluates the right-hand-side of a linear, Clohessy-
% Wiltshire propagator
%
% Author: Alessandro Cambielli
%
% Inputs:
%   state : relative state [6,1]
%   n : mean motion of mango [1,1]  
%
% Outputs:
%   dxdt : derivative of the state vector [6,1] 
%

% Initialize right-hand-side
dxdt = zeros(6,1);

x = state(1);
z = state(3);
xdot = state(4);
ydot = state(5);
zdot = state(6);

% position detivative is object's velocity
dxdt(1:3) = [xdot; ydot; zdot];   

% velocity derivative
dxdt(4:6) = [3*n^2*x+2*n*ydot; -2*n*xdot; -n^2*z];

end


