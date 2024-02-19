%% Spacecraft Guidance and Navigation - Periodic Orbit (2023/2024)
% Assignment:     2
% Exercise:       1 - Uncertainty propagation
% Author:         Alessandro Cambielli

% ID number:      10619322
%                 102615

%% REQUEST 1

% Clear memory workspace and set path
clearvars; close all; clc 

% Load SPICE kernels: 
cspice_furnsh('assignment02.tm');
format long g

%tsep = '2010-08-12-05:27:39.114 UTC'; 
tsep = cspice_str2et('2010-08-12-05:27:39.114 UTC'); % Ephemeris time [s]

mu = cspice_bodvrd('Earth','GM',1);  % Earth gravitational parameter [km^3/s^2]
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

r0_1 = [4622.232026629; 5399.3369588058; -0.0212138165769957];
v0_1 = [0.812221125483763; -0.721512914578826; 7.42665302729053];
x0_1 = [r0_1;v0_1];

r0_2 = [4621.69343340281; 5399.26386352847; -3.09039248714313];
v0_2 = [0.813960847513811; -0.719449862738607; 7.42706066911294];
x0_2 = [r0_2;v0_2];

P0 = [5.6e-7 3.5e-7 -7.1e-8 0 0 0;...
      3.5e-7 9.7e-7 7.6e-8 0 0 0;...
      -7.1e-8 7.6e-8 8.1e-8 0 0 0;...
      0 0 0 2.8e-11 0 0;...
       0 0 0 0 2.7e-11 0;...
       0 0 0 0 0 9.6e-12];

a1 = 1/(2/norm(r0_1)-norm(v0_1)^2/mu);
a2 = 1/(2/norm(r0_2)-norm(v0_2)^2/mu);
T = 2*pi*sqrt(a1^3/mu);
T2 = 2*pi*sqrt(a2^3/mu);

[~, xx1] = ode113(@(t,x) kepler(t,x,mu), [0 T], x0_1, options);
[~, xx2] = ode113(@(t,x) kepler(t,x,mu), [0 T2], x0_2, options);

tvec = tsep:T:tsep+10*T;

% Uncertainty propagation for satellite 1 and 2

% LinCov

mean_LC1 = zeros(length(tvec),6);
P_LC1 = zeros(6,6,length(tvec));
mean_LC1(1,:) = x0_1; 
P_LC1(:,:,1) = P0;

mean_LC2 = zeros(length(tvec),6);
P_LC2 = zeros(6,6,length(tvec));
mean_LC2(1,:) = x0_2; 
P_LC2(:,:,1) = P0;

for i = 2:length(tvec)
    [mean_LC1(i,:),P_LC1(:,:,i)] = LinCov(mean_LC1(i-1,:)',P_LC1(:,:,i-1),T);
    [mean_LC2(i,:),P_LC2(:,:,i)] = LinCov(mean_LC2(i-1,:)',P_LC2(:,:,i-1),T);
end

r_dispLC1 = zeros(length(mean_LC1),1);
v_dispLC1 = zeros(length(mean_LC1),1);
r_dispLC2 = zeros(length(mean_LC2),1);
v_dispLC2 = zeros(length(mean_LC2),1);

for i = 1:length(mean_LC1)
    r_dispLC1(i) = sqrt(trace(P_LC1(1:3,1:3,i)));
    v_dispLC1(i) = sqrt(trace(P_LC1(4:6,4:6,i)));
    r_dispLC2(i) = sqrt(trace(P_LC2(1:3,1:3,i)));
    v_dispLC2(i) = sqrt(trace(P_LC2(4:6,4:6,i)));
end

% UT

mean_UT1 = zeros(length(tvec),6);
P_UT1 = zeros(6,6,length(tvec));
mean_UT1(1,:) = x0_1; 
P_UT1(:,:,1) = P0;

mean_UT2 = zeros(length(tvec),6);
P_UT2 = zeros(6,6,length(tvec));
mean_UT2(1,:) = x0_2; 
P_UT2(:,:,1) = P0;

for i = 2:length(tvec)
    [mean_UT1(i,:),P_UT1(:,:,i)] = UT(mean_UT1(i-1,:)',P_UT1(:,:,i-1),T);
    [mean_UT2(i,:),P_UT2(:,:,i)] = UT(mean_UT2(i-1,:)',P_UT2(:,:,i-1),T);
end

r_dispUT1 = zeros(length(mean_UT1),1);
v_dispUT1 = zeros(length(mean_UT1),1);
r_dispUT2 = zeros(length(mean_UT2),1);
v_dispUT2 = zeros(length(mean_UT2),1);

for i = 1:length(mean_UT1)
    r_dispUT1(i) = sqrt(trace(P_UT1(1:3,1:3,i)));
    v_dispUT1(i) = sqrt(trace(P_UT1(4:6,4:6,i)));
    r_dispUT2(i) = sqrt(trace(P_UT2(1:3,1:3,i)));
    v_dispUT2(i) = sqrt(trace(P_UT2(4:6,4:6,i)));
end

%title('position standard deviation for satellite 1')
figure
hold on
grid on
plot(r_dispLC1,'ro')
plot(r_dispUT1,'b*')
xlim([0 12])
% Plot settings
set(gca,'FontSize',12)
legend('LinCov', 'UT', 'Location','northwest','Interpreter','latex','FontSize',14) 
xlabel('$period$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$position$ $standard$ $deviation$ [km]','Interpreter','latex','FontSize', 20)

%title('velocity standard deviation for satellite 1')
figure
hold on
grid on
plot(v_dispLC1,'ro')
plot(v_dispUT1,'b*')
xlim([0 12])
% Plot settings
set(gca,'FontSize',12)
legend('LinCov', 'UT', 'Location','northwest','Interpreter','latex','FontSize',14) 
xlabel('$period$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$velocity$ $standard$ $deviation$ [km/s]','Interpreter','latex','FontSize', 20)

%title('position standard deviation for satellite 2')
figure
hold on
grid on
plot(r_dispLC2,'ro')
plot(r_dispUT2,'b*')
xlim([0 12])
% Plot settings
set(gca,'FontSize',12)
legend('LinCov', 'UT', 'Location','northwest','Interpreter','latex','FontSize',14) 
xlabel('$period$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$position$ $standard$ $deviation$ [km]','Interpreter','latex','FontSize', 20)

%title('velocity standard deviation for satellite 2')
figure
hold on
grid on
plot(v_dispLC2,'ro')
plot(v_dispUT2,'b*')
xlim([0 12])
% Plot settings
set(gca,'FontSize',12)
legend('LinCov', 'UT', 'Location','northwest','Interpreter','latex','FontSize',14) 
xlabel('$period$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$velocity$ $standard$ $deviation$ [km/s]','Interpreter','latex','FontSize', 20)


%% REQUEST 2

deltar_LC = zeros(length(tvec),1);
deltar_UT = zeros(length(tvec),1);
Psum_LC = zeros(3,3,length(tvec));
Psum_UT = zeros(3,3,length(tvec));
limit_LC = zeros(length(tvec),1);
limit_UT = zeros(length(tvec),1);
counterLC = 0;
counterUT = 0;

for i = 1:length(tvec)
    deltar_LC(i) = norm(mean_LC1(i,1:3)-mean_LC2(i,1:3));
    Psum_LC(:,:,i) = P_LC1(1:3,1:3,i) + P_LC2(1:3,1:3,i);
    limit_LC(i) = 3*sqrt(eigs(Psum_LC(:,:,i),1));
    if deltar_LC(i) < limit_LC(i)
        break
    end
    counterLC = counterLC + 1;
end

for i = 1:length(tvec)
    deltar_UT(i) = norm(mean_UT1(i,1:3)-mean_UT2(i,1:3));
    Psum_UT(:,:,i) = P_UT1(1:3,1:3,i) + P_UT2(1:3,1:3,i);
    limit_UT(i) = 3*sqrt(eigs(Psum_UT(:,:,i),1));
    if deltar_UT(i) < limit_UT(i)
        break
    end
    counterUT = counterUT + 1;
end

% title('Revolution counter for LinCov')
figure
hold on
grid on
plot(deltar_LC(1:5),'r*')
plot(limit_LC(1:5),'b*')
xlim([0 6])
% Plot settings
set(gca,'FontSize',12)
legend('LinCov', 'limit', 'Location','northeast','Interpreter','latex','FontSize',14) 
xlabel('$period$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$relative$ $position$ $and$ $limit$ [km]','Interpreter','latex','FontSize', 20)

% title('Revolution counter for UT')
figure
hold on
grid on
plot(deltar_UT(1:5),'r*')
plot(limit_UT(1:5),'b*')
xlim([0 6])
% Plot settings
set(gca,'FontSize',12)
legend('UT', 'limit', 'Location','northeast','Interpreter','latex','FontSize',14) 
xlabel('$period$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$relative$ $position$ $and$ $limit$ [km]','Interpreter','latex','FontSize', 20)

%% REQUEST 3

n = 500;

mean_MCM1 = zeros(length(tvec),6);
P_MCM1 = zeros(6,6,length(tvec));
mean_MCM1(1,:) = x0_1; 
P_MCM1(:,:,1) = P0;

mean_MCM2 = zeros(length(tvec),6);
P_MCM2 = zeros(6,6,length(tvec));
mean_MCM2(1,:) = x0_2; 
P_MCM2(:,:,1) = P0;

for i = 2:length(tvec)
    [mean_MCM1(i,:),P_MCM1(:,:,i),ysamples1] = MCM(mean_MCM1(i-1,:)',P_MCM1(:,:,i-1),T,n);
    [mean_MCM2(i,:),P_MCM2(:,:,i),ysamples2] = MCM(mean_MCM2(i-1,:)',P_MCM2(:,:,i-1),T,n);
end

r_dispMCM1 = zeros(length(mean_MCM1),1);
v_dispMCM1 = zeros(length(mean_MCM1),1);
r_dispMCM2 = zeros(length(mean_MCM2),1);
v_dispMCM2 = zeros(length(mean_MCM2),1);

for i = 1:length(mean_MCM1)
    r_dispMCM1(i) = sqrt(trace(P_MCM1(1:3,1:3,i)));
    v_dispMCM1(i) = sqrt(trace(P_MCM1(4:6,4:6,i)));
    r_dispMCM2(i) = sqrt(trace(P_MCM2(1:3,1:3,i)));
    v_dispMCM2(i) = sqrt(trace(P_MCM2(4:6,4:6,i)));
end


% plots for LinCov, UT and MCM

eigr_LC1 = zeros(length(tvec),1);
eigr_LC2 = zeros(length(tvec),1);
eigv_LC1 = zeros(length(tvec),1);
eigv_LC2 = zeros(length(tvec),1);

eigr_UT1 = zeros(length(tvec),1);
eigr_UT2 = zeros(length(tvec),1);
eigv_UT1 = zeros(length(tvec),1);
eigv_UT2 = zeros(length(tvec),1);

eigr_MCM1 = zeros(length(tvec),1);
eigr_MCM2 = zeros(length(tvec),1);
eigv_MCM1 = zeros(length(tvec),1);
eigv_MCM2 = zeros(length(tvec),1);

for i = 1:length(tvec)

    eigr_LC1(i) = 3*sqrt(eigs(P_LC1(1:3,1:3,i),1));
    eigr_LC2(i) = 3*sqrt(eigs(P_LC2(1:3,1:3,i),1));
    eigv_LC1(i) = 3*sqrt(eigs(P_LC1(4:6,4:6,i),1));
    eigv_LC2(i) = 3*sqrt(eigs(P_LC1(4:6,4:6,i),1));

    eigr_UT1(i) = 3*sqrt(eigs(P_UT1(1:3,1:3,i),1));
    eigr_UT2(i) = 3*sqrt(eigs(P_UT2(1:3,1:3,i),1));
    eigv_UT1(i) = 3*sqrt(eigs(P_UT1(4:6,4:6,i),1));
    eigv_UT2(i) = 3*sqrt(eigs(P_UT1(4:6,4:6,i),1));

    eigr_MCM1(i) = 3*sqrt(eigs(P_MCM1(1:3,1:3,i),1));
    eigr_MCM2(i) = 3*sqrt(eigs(P_MCM2(1:3,1:3,i),1));
    eigv_MCM1(i) = 3*sqrt(eigs(P_MCM1(4:6,4:6,i),1));
    eigv_MCM2(i) = 3*sqrt(eigs(P_MCM1(4:6,4:6,i),1));
end

%title('position evolution for satellite 1')
figure
hold on
grid on
plot(eigr_LC1,'ro')
plot(eigr_UT1,'b*')
plot(eigr_MCM1,'m square')
xlim([0 12])
ylim([0 2.2])
% Plot settings
set(gca,'FontSize',12)
legend('LinCov', 'UT', 'Monte Carlo','Location','northwest',...
    'Interpreter','latex','FontSize',14) 
xlabel('$period$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$3 \sqrt{max(\lambda_1(P_{r,1})}$ [km]','Interpreter','latex','FontSize', 20)


%title('velocity evolution for satellite 1')
figure
hold on
grid on
plot(eigv_LC1,'ro')
plot(eigv_UT1,'b*')
plot(eigv_MCM1,'m square')
xlim([0 12])
ylim([0 2.3e-3])
% Plot settings
set(gca,'FontSize',12)
legend('LinCov', 'UT', 'Monte Carlo','Location','northwest',...
    'Interpreter','latex','FontSize',14) 
xlabel('$period$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$3 \sqrt{max(\lambda_1(P_{v,1})}$ [km/s]','Interpreter','latex','FontSize', 20)


%title('position evolution for satellite 2')
figure
hold on
grid on
plot(eigr_LC2,'ro')
plot(eigr_UT2,'b*')
plot(eigr_MCM2,'m square')
xlim([0 12])
ylim([0 2.2])
% Plot settings
set(gca,'FontSize',12)
legend('LinCov', 'UT', 'Monte Carlo','Location','northwest',...
    'Interpreter','latex','FontSize',14) 
xlabel('$period$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$3 \sqrt{max(\lambda_2(P_{r,2})}$ [km]','Interpreter','latex','FontSize', 20)


%title('velocity evolution for satellite 2')
figure
hold on
grid on
plot(eigv_LC2,'ro')
plot(eigv_UT2,'b*')
plot(eigv_MCM2,'m square')
xlim([0 12])
ylim([0 2.3e-3])
% Plot settings
set(gca,'FontSize',12)
legend('LinCov', 'UT', 'Monte Carlo','Location','northwest',...
    'Interpreter','latex','FontSize',14) 
xlabel('$period$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$3 \sqrt{max(\lambda_2(P_{v,2})}$ [km/s]','Interpreter','latex','FontSize', 20)

% PLOTTING THE ELLIPSES FOR SAT1

% rotation to LVLH frame
i = r0_1/norm(r0_1);
k = cross(r0_1,v0_1)/norm(cross(r0_1,v0_1));
j = cross(k,i);
R = [i, j, k]';

ysamples_rot1 = zeros(n,3);

for i = 1:n
    ysamples_rot1(i,:) = R*ysamples1(i,1:3)';
end

mean_LC1_rot1 = R*mean_LC1(1,1:3)';
mean_UT1_rot1 = R*mean_UT1(1,1:3)';
mean_MCM1_rot1 = R*mean_MCM1(1,1:3)';

mean_LC1_rot2 = R*mean_LC1(end,1:3)';
mean_UT1_rot2 = R*mean_UT1(end,1:3)';
mean_MCM1_rot2 = R*mean_MCM1(end,1:3)';

P_LC1_rot1 = R*P_LC1(1:3,1:3,1)*R';
P_UT1_rot1 = R*P_UT1(1:3,1:3,1)*R';
P_MCM1_rot1 = R*P_MCM1(1:3,1:3,1)*R';

P_LC1_rot2 = R*P_LC1(1:3,1:3,end)*R';
P_UT1_rot2 = R*P_UT1(1:3,1:3,end)*R';
P_MCM1_rot2 = R*P_MCM1(1:3,1:3,end)*R';


% title('error ellipses for LinCov')
figure
hold on
grid on
plot(ysamples_rot1(:,1),ysamples_rot1(:,2),'.','Color',[.7 .7 .7])
ellipse(mean_LC1_rot1(1:2),P_LC1_rot1(1:2,1:2));
ellipse(mean_LC1_rot2(1:2),P_LC1_rot2(1:2,1:2));
scatter(mean_LC1_rot1(1),mean_LC1_rot1(2),80,...
    'MarkerFaceColor','none')
scatter(mean_LC1_rot2(1),mean_LC1_rot2(2),80,...
    'MarkerFaceColor','none')
% Plot settings
set(gca,'FontSize',12)
legend('MC samples','first revolution ellipse','final revolution ellipse',...
    'first revolution mean','final revolution mean','Location','southwest',...
    'Interpreter','latex','FontSize',14) 
xlabel('$i$ [km]','Interpreter','latex','FontSize', 20)
ylabel('$j$ [km]','Interpreter','latex','FontSize', 20)


% title('error ellipses for UT')
figure
hold on
grid on
plot(ysamples_rot1(:,1),ysamples_rot1(:,2),'.','Color',[.7 .7 .7])
ellipse(mean_UT1_rot1(1:2),P_UT1_rot1(1:2,1:2));
ellipse(mean_UT1_rot2(1:2),P_UT1_rot2(1:2,1:2));
scatter(mean_UT1_rot1(1),mean_UT1_rot1(2),80,...
    'MarkerFaceColor','none')
scatter(mean_UT1_rot2(1),mean_UT1_rot2(2),80,...
    'MarkerFaceColor','none')
% Plot settings
set(gca,'FontSize',12)
legend('MC samples','first revolution ellipse','final revolution ellipse',...
    'first revolution mean','final revolution mean','Location','southwest',...
    'Interpreter','latex','FontSize',14) 
xlabel('$i$ [km]','Interpreter','latex','FontSize', 20)
ylabel('$j$ [km]','Interpreter','latex','FontSize', 20)


% title('error ellipses for MCM')
figure
hold on
grid on
plot(ysamples_rot1(:,1),ysamples_rot1(:,2),'.','Color',[.7 .7 .7])
ellipse(mean_MCM1_rot1(1:2),P_MCM1_rot1(1:2,1:2));
ellipse(mean_MCM1_rot2(1:2),P_MCM1_rot2(1:2,1:2));
scatter(mean_MCM1_rot1(1),mean_MCM1_rot1(2),80,...
    'MarkerFaceColor','none')
scatter(mean_MCM1_rot2(1),mean_MCM1_rot2(2),80,...
    'MarkerFaceColor','none')
% Plot settings
set(gca,'FontSize',12)
legend('MC samples','first revolution ellipse','final revolution ellipse',...
    'first revolution mean','final revolution mean','Location','southwest',...
    'Interpreter','latex','FontSize',14) 
xlabel('$i$ [km]','Interpreter','latex','FontSize', 20)
ylabel('$j$ [km]','Interpreter','latex','FontSize', 20)

% PLOTTING THE ELLIPSES FOR SAT2

% rotation to LVLH frame
i = r0_2/norm(r0_2);
k = cross(r0_2,v0_2)/norm(cross(r0_2,v0_2));
j = cross(k,i);
R = [i, j, k]';

ysamples_rot2 = zeros(n,3);

for i = 1:n
    ysamples_rot2(i,:) = R*ysamples2(i,1:3)';
end

mean_LC1_rot1 = R*mean_LC2(1,1:3)';
mean_UT1_rot1 = R*mean_UT2(1,1:3)';
mean_MCM1_rot1 = R*mean_MCM2(1,1:3)';

mean_LC1_rot2 = R*mean_LC2(end,1:3)';
mean_UT1_rot2 = R*mean_UT2(end,1:3)';
mean_MCM1_rot2 = R*mean_MCM2(end,1:3)';

P_LC1_rot1 = R*P_LC2(1:3,1:3,1)*R';
P_UT1_rot1 = R*P_UT2(1:3,1:3,1)*R';
P_MCM1_rot1 = R*P_MCM2(1:3,1:3,1)*R';

P_LC1_rot2 = R*P_LC2(1:3,1:3,end)*R';
P_UT1_rot2 = R*P_UT2(1:3,1:3,end)*R';
P_MCM1_rot2 = R*P_MCM2(1:3,1:3,end)*R';


% title('error ellipses for LinCov')
figure
hold on
grid on
plot(ysamples_rot2(:,1),ysamples_rot2(:,2),'.','Color',[.7 .7 .7])
ellipse(mean_LC1_rot1(1:2),P_LC1_rot1(1:2,1:2));
ellipse(mean_LC1_rot2(1:2),P_LC1_rot2(1:2,1:2));
scatter(mean_LC1_rot1(1),mean_LC1_rot1(2),80,...
    'MarkerFaceColor','none')
scatter(mean_LC1_rot2(1),mean_LC1_rot2(2),80,...
    'MarkerFaceColor','none')
% Plot settings
set(gca,'FontSize',12)
legend('MC samples','first revolution ellipse','final revolution ellipse',...
    'first revolution mean','final revolution mean','Location','southwest',...
    'Interpreter','latex','FontSize',14) 
xlabel('$i$ [km]','Interpreter','latex','FontSize', 20)
ylabel('$j$ [km]','Interpreter','latex','FontSize', 20)


% title('error ellipses for UT')
figure
hold on
grid on
plot(ysamples_rot2(:,1),ysamples_rot2(:,2),'.','Color',[.7 .7 .7])
ellipse(mean_UT1_rot1(1:2),P_UT1_rot1(1:2,1:2));
ellipse(mean_UT1_rot2(1:2),P_UT1_rot2(1:2,1:2));
scatter(mean_UT1_rot1(1),mean_UT1_rot1(2),80,...
    'MarkerFaceColor','none')
scatter(mean_UT1_rot2(1),mean_UT1_rot2(2),80,...
    'MarkerFaceColor','none')
% Plot settings
set(gca,'FontSize',12)
legend('MC samples','first revolution ellipse','final revolution ellipse',...
    'first revolution mean','final revolution mean','Location','southwest',...
    'Interpreter','latex','FontSize',14) 
xlabel('$i$ [km]','Interpreter','latex','FontSize', 20)
ylabel('$j$ [km]','Interpreter','latex','FontSize', 20)


% title('error ellipses for MCM')
figure
hold on
grid on
plot(ysamples_rot2(:,1),ysamples_rot2(:,2),'.','Color',[.7 .7 .7])
ellipse(mean_MCM1_rot1(1:2),P_MCM1_rot1(1:2,1:2));
ellipse(mean_MCM1_rot2(1:2),P_MCM1_rot2(1:2,1:2));
scatter(mean_MCM1_rot1(1),mean_MCM1_rot1(2),80,...
    'MarkerFaceColor','none')
scatter(mean_MCM1_rot2(1),mean_MCM1_rot2(2),80,...
    'MarkerFaceColor','none')
% Plot settings
set(gca,'FontSize',12)
legend('MC samples','first revolution ellipse','final revolution ellipse',...
    'first revolution mean','final revolution mean','Location','southwest',...
    'Interpreter','latex','FontSize',14) 
xlabel('$i$ [km]','Interpreter','latex','FontSize', 20)
ylabel('$j$ [km]','Interpreter','latex','FontSize', 20)



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

function dstate = kepler_STM(~, state, mu)
% The function provides the vectorial field for the circular restricted
% three-body problem merging the state with the computation of the state
% transition matrix
%
% Author: Alessandro Cambielli
%  
% Inputs:
%   state : cartesian state and STM [42,1]
%   mu : gravitational parameter [1,1]    
%
% Outputs:
%   dstate : derivative of cartesian state and STM [42,1] 
%

%Initialize
dstate = zeros(42,1);

x = state(1);
y = state(2);
z = state(3);
rr = [x; y; z]; % column
STM = reshape(state(7:42),6,6)';   %From equations to STM

% Compute distance
dist = sqrt(x^2 + y^2 + z^2);

% Compute the gravitational acceleration using Newton's law
aa_grav = - mu * rr /(dist^3);

% Compute the derivative of the flow
dfdv = 3*mu/dist^5*(rr*rr') - mu/dist^3*eye(3);

% Assemble the matrix A(t) = dfdx
A = [zeros(3), eye(3);...
           dfdv, zeros(3)];

% STM derivative
dSTM = A*STM;

% assemble derivatives of cartesian state and STM
dstate(1:3) = state(4:6);
dstate(4:6) = aa_grav;
dstate(7:12) = dSTM(1,1:6)';
dstate(13:18) = dSTM(2,1:6)';
dstate(19:24) = dSTM(3,1:6)';
dstate(25:30) = dSTM(4,1:6)';
dstate(31:36) = dSTM(5,1:6)';
dstate(37:42) = dSTM(6,1:6)';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mean,P] = LinCov(x0,P0,T)
% The function simulates a linearized approach to propagate initial mean
% state and covariance.
%
% Author: Alessandro Cambielli
%
% Inputs:
%   x0 : initial mean state [6,1] 
%   P0 : initial covariance matrix [6,6] 
%   T : propagation time [1,1] 
%
% Outputs:
%   mean : final mean state [6,1] 
%   P : final covariance matrix [6,6]
%

mu = cspice_bodvrd('Earth','GM',1);
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

[~, xx] = ode113(@(t,x) kepler_STM(t,x,mu), [0 T],...
          [x0; reshape(eye(6),[],1)], options);
mean = xx(end,1:6);
STM = (reshape(xx(end,7:42),6,6))';
P = STM*P0*STM';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mean,P] = UT(mean0,P0,T)
% The function propagates initial mean state and covariance exploiting the
% unscented transform.
%
% Author: Alessandro Cambielli
%
% Inputs:
%   mean0 : initial mean state [6,1] 
%   P0 : initial covariance matrix [6,6] 
%   T : propagation time [1,1] 
%
% Outputs:
%   mean : final mean state [6,1] 
%   P : final covariance matrix [6,6]
%

mu = cspice_bodvrd('Earth','GM',1);
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

n = 6;
alpha = 0.1;
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
    [~, xx] = ode113(@(t,x) kepler(t,x,mu), [0 T], chi(:,i), options);
    jota(:,i) = xx(end,:)';
    mean = mean + w_m(i) * jota(:,i);
end

for k = 1:(2*n+1)
    P = P + w_c(k) * (jota(:,k)-mean) * (jota(:,k)-mean)';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mean_MCM,P_MCM,ysamples] = MCM(x0,P0,T,n)
% The function performs a Monte Carlo analysis to propagate initial mean 
% state and covariance.
%
% Author: Alessandro Cambielli
%
% Inputs:
%   x0 : initial mean state [6,1] 
%   P0 : initial covariance matrix [6,6] 
%   T : propagation time [1,1]
%   n : number of samples drawn from the initial covariance [1,1]
%
% Outputs:
%   mean_MCM : final mean state [6,1] 
%   P_MCM : final covariance matrix [6,6]
%   ysamples : propagated samples [n,6]
%

mu = cspice_bodvrd('Earth','GM',1);
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);

rng('default')  % For reproducibility

ysamples = zeros(n,6);

xsamples = mvnrnd(x0,P0,n);

for i = 1:length(xsamples)
    [~, prop] = ode113(@(t,x) kepler(t,x,mu),[0 T],...
        xsamples(i,:)', options);
    ysamples(i,:) = prop(end,:);
end

mean_MCM = mean(ysamples);
P_MCM = cov(ysamples);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = ellipse(mean,covariance)
% The function draws the 3 sigma error ellipse derived from the mean
% position spread.
%
% Author: Alessandro Cambielli
%
% Inputs:
%   mean : mean state [3,1] 
%   covariance : covariance matrix [3,3] 
%
% Outputs:
%   value : expression describing the event [1,1]
%   isterminal : 1 if the integration is to terminate when the event occurs
%   direction : -1 locates only zeros where the event function is decreasing
%

points = 200; % number of plotted points on the ellipse

alpha  = 2*pi/points*(0:points);
circle = [cos(alpha);sin(alpha)];
sigma = 3; 
P = covariance;
[R,D] = svd(P); % single value decomposition
d = sqrt(D);
ellip = sigma*R*d*circle;

% retrieve points
x = mean(1)+ellip(1,:);
y = mean(2)+ellip(2,:);

% ellipse plot
plot(x,y,'LineWidth',2);

end

