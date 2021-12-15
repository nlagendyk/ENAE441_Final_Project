%% ENAE441_Group Project
% Group 8, Section #0101
% Team Members: Ronak Chawla, Prashanta Aryal, Nico Lagendyk, Ryan Quigley
clear;
clc;
close all;

% Please reference the several functions we created and used at the bottom of this .m file
%% Part 2: Data
% Load in the data:
load('opt2satBset4');
load('opt3satBset4');

%% Part 3.1: Initial Orbit Determination

mkdir('ResidualGraphs_IOD');

% IOD: Early Night 1 Points
night1_early.timespan = [10:4:110]; % Observation data set used
night1_early.vector_set = [10 15 20;  10 55 110;  100 105 110];

[night1_early.IOD.early, night2_early.IOD_prop.early] = IOD(opt2satBset4,opt3satBset4,night1_early.timespan,night1_early.vector_set(1,:),'night1 early IOD early');
[night1_early.IOD.spread, night2_early.IOD_prop.spread] = IOD(opt2satBset4,opt3satBset4,night1_early.timespan,night1_early.vector_set(2,:),'night1 early IOD spread');
[night1_early.IOD.late, night2_early.IOD_prop.late] = IOD(opt2satBset4,opt3satBset4,night1_early.timespan,night1_early.vector_set(3,:),'night1 early IOD late');

% IOD: Mid Night 1 Points
night1_mid.timespan = [350:4:450];
night1_mid.vector_set = [350 355 360;  350 400 450;  440 445 450];
                     
[night1_mid.IOD.early, night2_mid.IOD_prop.early] = IOD(opt2satBset4,opt3satBset4,night1_mid.timespan,night1_mid.vector_set(1,:),'night1 mid IOD early');
[night1_mid.IOD.spread, night2_mid.IOD_prop.spread] = IOD(opt2satBset4,opt3satBset4,night1_mid.timespan,night1_mid.vector_set(2,:),'night1 mid IOD spread');
[night1_mid.IOD.late, night2_mid.IOD_prop.late] = IOD(opt2satBset4,opt3satBset4,night1_mid.timespan,night1_mid.vector_set(3,:),'night1 mid IOD late');

% IOD: Late Night 1Points
night1_late.timespan = [600:4:700];
night1_late.vector_set = [600 605 610;  600 650 699;  690 695 700];
                     
[night1_late.IOD.early, night2_late.IOD_prop.early] = IOD(opt2satBset4,opt3satBset4,night1_late.timespan,night1_late.vector_set(1,:),'night1 late IOD early');
[night1_late.IOD.spread, night2_late.IOD_prop.spread]  = IOD(opt2satBset4,opt3satBset4,night1_late.timespan,night1_late.vector_set(2,:),'night1 late IOD spread');
[night1_late.IOD.late, night2_late.IOD_prop.late] = IOD(opt2satBset4,opt3satBset4,night1_late.timespan,night1_late.vector_set(3,:),'night1 late IOD late');

% IOD: Night 1 Spread out Points
night1_spread.timespan = [10:10:700];
night1_spread.vector_set = [10 50 100;  10 350 700;  600 650 700];
                     
[night1_spread.IOD.early, night2_spread.IOD_prop.early] = IOD(opt2satBset4,opt3satBset4,night1_spread.timespan,night1_spread.vector_set(1,:),'night1 spread IOD early');
[night1_spread.IOD.spread, night2_spread.IOD_prop.spread]  = IOD(opt2satBset4,opt3satBset4,night1_spread.timespan,night1_spread.vector_set(2,:),'night1 spread IOD spread');
[night1_spread.IOD.late, night2_spread.IOD_prop.late] = IOD(opt2satBset4,opt3satBset4,night1_spread.timespan,night1_spread.vector_set(3,:),'night1 spread IOD late');

% Compare the RMS values and find the best initial estimate from IOD:
observation_set = ["night1_early.IOD.early"; "night1_early.IOD.spread"; "night1_early.IOD.late";
                   "night1_mid.IOD.early"; "night1_mid.IOD.spread"; "night1_mid.IOD.late";
                   "night1_late.IOD.early"; "night1_late.IOD.spread"; "night1_late.IOD.late";
                   "night1_spread.IOD.early"; "night1_spread.IOD.spread"; "night1_spread.IOD.late"];

best_rms_value = [night1_early.IOD.early.best.rms; night1_early.IOD.spread.best.rms; night1_early.IOD.late.best.rms;
                  night1_mid.IOD.early.best.rms; night1_mid.IOD.spread.best.rms; night1_mid.IOD.late.best.rms;
                  night1_late.IOD.early.best.rms; night1_late.IOD.spread.best.rms; night1_late.IOD.late.best.rms;
                  night1_spread.IOD.early.best.rms; night1_spread.IOD.spread.best.rms; night1_spread.IOD.late.best.rms];

IOD_compare = table(observation_set,best_rms_value)

[IOD_bestset,j] = min(IOD_compare.best_rms_value);

fprintf('The best rms value we get from the inital orbit determination is %0.4f\n', IOD_bestset)
fprintf('which is the best rms value from the %s observation set\n',IOD_compare.observation_set(j))

save('IOD_workspace.mat') 

%% Part 3.2: Estimation
% clear;
% clc;
% load('IOD_workspace.mat');

force_model_4x4_1m2_cr1 = force_model_third_body(4, 4, 0, 0, 1, 1, 1000);
oapchile = make_station("OAP-Chile",-30.1428030000, -70.6945280000, 1500.000000);
fourcols = ["observation_number" "datetime" "azimuth_deg" "elevation_deg"];
observer_lla = latlonalt_deg(opt2satBset4.site_latitude_deg(1), opt2satBset4.site_longitude_deg(1), opt2satBset4.site_altitude_m(1));

% The best initial estimate is:
initial_est = night1_spread.IOD.late.t1.initial_est;
fprintf('\nThe best initial estimate is night1_spread.IOD.late.t1.initial_est\n\n')

%% Part 3.2.2: Study

% Night 1 Observation sets for Orbit determination:
night1_early.pts = opt2satBset4(night1_early.timespan,fourcols);
night1_mid.pts = opt2satBset4(night1_mid.timespan,fourcols);
night1_late.pts = opt2satBset4(night1_late.timespan,fourcols);
night1_spread.pts = opt2satBset4(night1_spread.timespan,fourcols);

%% Study 1. Vary the observation data set, holding the propagator constant:

% Orbit Determination:
od_night1_early = determine_orbit(initial_est, oapchile, night1_early.pts, force_model_4x4_1m2_cr1);
od_night1_mid = determine_orbit(initial_est, oapchile, night1_mid.pts, force_model_4x4_1m2_cr1);
od_night1_late = determine_orbit(initial_est, oapchile, night1_late.pts, force_model_4x4_1m2_cr1);
od_night1_spread = determine_orbit(initial_est, oapchile, night1_spread.pts, force_model_4x4_1m2_cr1);

% Night 2 Observed values and times to Propogate to:
night2_early.pts = opt3satBset4([10:4:110],fourcols);
night2_mid.pts = opt3satBset4([350:4:450],fourcols);
night2_late.pts = opt3satBset4([600:4:700],fourcols);
night2_spread.pts = opt3satBset4([10:10:700],fourcols);

% Propogating while varying observation sets and holding force model constant:
night2_early.prop = propagate_to_times(od_night1_early.estimated, night2_early.pts.datetime, force_model_4x4_1m2_cr1);
night2_mid.prop = propagate_to_times(od_night1_mid.estimated, night2_mid.pts.datetime, force_model_4x4_1m2_cr1);
night2_late.prop = propagate_to_times(od_night1_late.estimated, night2_late.pts.datetime, force_model_4x4_1m2_cr1);
night2_spread.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_4x4_1m2_cr1);

% Converting propogated position and velocities into azimuths and
% elevations to compare these predicted values with the observed:
night2_early.prop = ADD_aer(night2_early.prop, observer_lla);
night2_mid.prop = ADD_aer(night2_mid.prop, observer_lla);
night2_late.prop = ADD_aer(night2_late.prop, observer_lla);
night2_spread.prop = ADD_aer(night2_spread.prop, observer_lla);

% Calculate RMS: 

night2_early.pts.epoch = night2_early.pts.datetime;
night2_mid.pts.epoch = night2_mid.pts.datetime;
night2_late.pts.epoch = night2_late.pts.datetime;
night2_spread.pts.epoch = night2_spread.pts.datetime;

% Plot Residuals of Night 1 Observation Sets Propagated to Night 2 vs. Observed Night 2 Data:

subplot(2,2,1)
hold on
night1_early.prop.RMS = RMS(night2_early.pts, night2_early.prop);
title('Early Observation Set')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

subplot(2,2,2)
hold on
night1_mid.prop.RMS = RMS(night2_mid.pts, night2_mid.prop);
title('Middle Observation Set')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

subplot(2,2,3)
hold on
night1_late.prop.RMS = RMS(night2_late.pts, night2_late.prop);
title('Late Observation Set')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

subplot(2,2,4)
hold on
night1_spread.prop.RMS = RMS(night2_spread.pts, night2_spread.prop);
title('Spread Observation Set')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

sgtitle({'Varying the Observation Sets:', 'Residuals of Night 1 Obs Sets Propagated to Night 2 vs. Obs Night 2 Data'})
sgt.FontSize = 24;

mkdir('Study_graphs');

% Save the figures to local directory:
if ismac
    saveas(gcf,[pwd sprintf('/Study_graphs/%s','Changing the Observation Sets')]) % for mac
elseif ispc
    saveas(gcf,[pwd sprintf('\\Study_graphs\\%s','Changing the Observation Sets')]) % for win
end


% Compare the RMS values and find the best Observation Set from Night 1:
observation_set = ["night1_early";
                   "night1_mid";
                   "night1_late";
                   "night1_spread"];

best_rms_value = [night1_early.prop.RMS.value;
                  night1_mid.prop.RMS.value;
                  night1_late.prop.RMS.value;
                  night1_spread.prop.RMS.value];

obs_compare = table(observation_set,best_rms_value)

[obs_best_set,k] = min(obs_compare.best_rms_value);

fprintf('The best rms value we get from varying the observation data sets is %0.8f\n', obs_best_set)
fprintf('which is from the %s observation set\n',obs_compare.observation_set(k))

save('study1_workspace.mat')

%% Study 2: Vary the force model, holding the observation data set constant:

%% Study 2.1: Changing the gravity force model:
% Gravity forces to test include two body, 2 × 0, 2 × 2, and 20 × 20.  

% Force models defined:
force_model_4x4 = force_model(4, 4, 0, 0, 1, 1, 6000);
force_model_2x0 = force_model(2, 0, 0, 0, 1, 1, 6000);
force_model_2x2 = force_model(2, 2, 0, 0, 1, 1, 6000);
force_model_20x20 = force_model(20, 20, 0, 0, 1, 1, 6000);  

% Propogating while varying the force models and holding the observation set constant:
% For this study, we chose to use the od_night1_spread observation set.
night2_4x4.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_4x4);
night2_2x0.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_2x0);
night2_2x2.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_2x2);
night2_20x20.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_20x20);

% Converting propogated position and velocities into azimuths and
% elevations to compare these predicted values with the observed:
night2_4x4.prop = ADD_aer(night2_4x4.prop, observer_lla);
night2_2x0.prop = ADD_aer(night2_2x0.prop, observer_lla);
night2_2x2.prop = ADD_aer(night2_2x2.prop, observer_lla);
night2_20x20.prop = ADD_aer(night2_20x20.prop, observer_lla);

% Calculate RMS: 

% Plot Residuals of Night 1 Observation Sets Propagated to Night 2 vs. Observed Night 2 Data:

subplot(2,2,1)
hold on
night2_4x4.RMS = RMS(night2_spread.pts, night2_4x4.prop);
title('night2 4x4 force model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

subplot(2,2,2)
hold on
night2_2x0.RMS = RMS(night2_spread.pts, night2_2x0.prop);
title('night2 2x0 force model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

subplot(2,2,3)
hold on
night2_2x2.RMS = RMS(night2_spread.pts, night2_2x2.prop);
title('night2 2x2 force model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

subplot(2,2,4)
hold on
night2_20x20.RMS = RMS(night2_spread.pts, night2_20x20.prop);
title('night2 20x20 force model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

sgtitle({'Changing the Gravity Force Model:', 'Residuals of Night 1 Obs Sets Propagated to Night 2 vs. Obs Night 2 Data'})
sgt.FontSize = 24;

% Save the figures to local directory:
if ismac
    saveas(gcf,[pwd sprintf('/Study_graphs/%s','Changing the Gravity Force Model')]) % for mac
elseif ispc
    saveas(gcf,[pwd sprintf('\\Study_graphs\\%s','Changing the Gravity Force Model')]) % for win
end

% Compare the RMS values and find the best force model:
observation_set = ["night2_4x4";
                   "night2_2x0";
                   "night2_2x2";
                   "night2_20x20"];

best_rms_value = [night2_4x4.RMS.value;
                  night2_2x0.RMS.value;
                  night2_2x2.RMS.value;
                  night2_20x20.RMS.value];

fm_compare = table(observation_set,best_rms_value)

[fm_best_set,l] = min(fm_compare.best_rms_value);

fprintf('The best rms value we get from varying the two-body force models is %0.8f\n', fm_best_set)
fprintf('which is from the %s observation set\n',fm_compare.observation_set(l))

save('study2_1_workspace.mat')

%% Study 2.2: Changing solar pressure parameters:

%% Varying Mass
format long
force_model_mass1 = force_model_third_body(2, 0, 0, 0, 1, 1.2, 1000); 
force_model_mass2 = force_model_third_body(2, 0, 0, 0, 1, 1.2, 5000); 
force_model_mass3 = force_model_third_body(2, 0, 0, 0, 1, 1.2, 9000);

% Propogating while varying the force models and holding the observation set constant:
% For this study, we chose to use the od_night1_spread observation set.
night2_mass1.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_mass1);
night2_mass2.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_mass2);
night2_mass3.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_mass3);

% Converting propogated position and velocities into azimuths and
% elevations to compare these predicted values with the observed:
night2_mass1.prop = ADD_aer(night2_mass1.prop, observer_lla);
night2_mass2.prop = ADD_aer(night2_mass2.prop, observer_lla);
night2_mass3.prop = ADD_aer(night2_mass3.prop, observer_lla);

% Calculate RMS: 

% Plot Residuals of Night 1 Observation Sets Propagated to Night 2 vs. Observed Night 2 Data:
subplot(3,1,1)
hold on
night2_mass1.RMS = RMS(night2_spread.pts, night2_mass1.prop);
title('night2 mass1 = 1000 kgforce model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

subplot(3,1,2)
hold on
night2_mass2.RMS = RMS(night2_spread.pts, night2_mass2.prop);
title('night2 mass2 = 5000 kg force model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

subplot(3,1,3)
hold on
night2_mass3.RMS = RMS(night2_spread.pts, night2_mass3.prop);
title('night2 mass3 = 9000 kg force model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

sgtitle({'Varying Mass:',' Residuals of Night 1 Obs Sets Propagated to Night 2 vs. Obs Night 2 Data'})
sgt.FontSize = 24;

% Save the figures to local directory:
if ismac
    saveas(gcf,[pwd sprintf('/Study_graphs/%s','Changing Solar Pressure Parameters Varying Mass')]) % for mac
elseif ispc
    saveas(gcf,[pwd sprintf('\\Study_graphs\\%s','Changing Solar Pressure Parameters Varying Mass')]) % for win
end

% Compare the RMS values and find the best force model:
observation_set = ["night2_mass1: 1000 kg";
                   "night2_mass2: 5000 kg";
                   "night2_mass3: 9000 kg"];

best_rms_value = [night2_mass1.RMS.value;
                  night2_mass2.RMS.value;
                  night2_mass3.RMS.value];

mass_compare = table(observation_set,best_rms_value)

[m_best_set,m] = min(mass_compare.best_rms_value);

fprintf('The best rms value we get from varying the mass parameter is %0.8f\n', m_best_set)
fprintf('which is from the %s observation set\n',mass_compare.observation_set(m))

%% Varying Surface area

%Force model varying surface area
force_model_Area1 = force_model_third_body(2, 0, 0, 0, .5, 1.2, 6000); 
force_model_Area2 = force_model_third_body(2, 0, 0, 0, 2, 1.2, 6000); 
force_model_Area3 = force_model_third_body(2, 0, 0, 0, 4, 1.2, 6000); 

% Propogating while varying the force models and holding the observation set constant:
% For this study, we chose to use the od_night1_spread observation set.
night2_Area1.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_Area1);
night2_Area2.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_Area2);
night2_Area3.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_Area3);

% Converting propogated position and velocities into azimuths and
% elevations to compare these predicted values with the observed:
night2_Area1.prop = ADD_aer(night2_Area1.prop, observer_lla);
night2_Area2.prop = ADD_aer(night2_Area2.prop, observer_lla);
night2_Area3.prop = ADD_aer(night2_Area3.prop, observer_lla);

% Calculate RMS: 

% Plot Residuals of Night 1 Observation Sets Propagated to Night 2 vs. Observed Night 2 Data:
figure
subplot(3,1,1)
hold on
night2_Area1.RMS = RMS(night2_spread.pts, night2_Area1.prop);
title('night2 Area1 = 1.5 m^2 force model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

subplot(3,1,2)
hold on
night2_Area2.RMS = RMS(night2_spread.pts, night2_Area2.prop);
title('night2 Area2 = 2 m^2 force model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

subplot(3,1,3)
hold on
night2_Area3.RMS = RMS(night2_spread.pts, night2_Area3.prop);
title('night2 Area3 = 4 m^2 force model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

sgtitle({'Varying Surface Area:', 'Residuals of Night 1 Obs Sets Propagated to Night 2 vs. Obs Night 2 Data'})
sgt.FontSize = 24;

% Save the figures to local directory:
if ismac
    saveas(gcf,[pwd sprintf('/Study_graphs/%s','Changing Solar Pressure Parameters Varying Surface Area')]) % for mac
elseif ispc
    saveas(gcf,[pwd sprintf('\\Study_graphs\\%s','Changing Solar Pressure Parameters Varying Surface Area')]) % for win
end

% Compare the RMS values and find the best force model:
observation_set = ["night2_Area1: 0.5 m^2";
                   "night2_Area2: 2 m^2";
                   "night2_Area3: 4 m^2"];

best_rms_value = [night2_Area1.RMS.value;
                  night2_Area2.RMS.value;
                  night2_Area3.RMS.value];

area_compare = table(observation_set,best_rms_value)

[area_best_set,n] = min(area_compare.best_rms_value);

fprintf('The best rms value we get from varying the cross-sectional area parameter is %0.8f\n', area_best_set)
fprintf('which is from the %s observation set\n',area_compare.observation_set(n))

%% Varying Coefficient of Solar Reflectivity

%Force model varying coefficient of solar reflectivity
force_model_cr1 = force_model_third_body(2, 0, 0, 0, 1, 0.6, 6000); 
force_model_cr2 = force_model_third_body(2, 0, 0, 0, 1, 1.2, 6000); 
force_model_cr3 = force_model_third_body(2, 0, 0, 0, 1, 1.8, 6000); 

% Propogating while varying the force models and holding the observation set constant:
% For this study, we chose to use the od_night1_spread observation set.
night2_cr1.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_cr1);
night2_cr2.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_cr2);
night2_cr3.prop = propagate_to_times(od_night1_spread.estimated, night2_spread.pts.datetime, force_model_cr3);

% Converting propogated position and velocities into azimuths and
% elevations to compare these predicted values with the observed:
night2_cr1.prop = ADD_aer(night2_cr1.prop, observer_lla);
night2_cr2.prop = ADD_aer(night2_cr2.prop, observer_lla);
night2_cr3.prop = ADD_aer(night2_cr3.prop, observer_lla);

%% Calculate RMS: 

% Plot Residuals of Night 1 Observation Sets Propagated to Night 2 vs. Observed Night 2 Data:
figure
subplot(3,1,1)
hold on
night2_cr1.RMS = RMS(night2_spread.pts, night2_cr1.prop);
title('night2 reflectivity1 = 0.6 force model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

subplot(3,1,2)
hold on
night2_cr2.RMS = RMS(night2_spread.pts, night2_cr2.prop);
title('night2 reflectivity2 = 1.2 force model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

subplot(3,1,3)
hold on
night2_cr3.RMS = RMS(night2_spread.pts, night2_cr3.prop);
title('night2 reflectivity3 = 1.8 force model')
xlabel('Date Time')
ylabel('Residual (deg)')
legend('Azimuth','Elevation')

sgtitle({'Varying Coefficient of Reflectivity:','Residuals of Night 1 Obs Sets Propagated to Night 2 vs. Obs Night 2 Data'})
sgt.FontSize = 24;

% Save the figures to local directory:
if ismac
    saveas(gcf,[pwd sprintf('/Study_graphs/%s','Changing Solar Pressure Parameters Varying Coefficient of Reflectivity')]) % for mac
elseif ispc
    saveas(gcf,[pwd sprintf('\\Study_graphs\\%s','Changing Solar Pressure Parameters Varying Coefficient of Reflectivity')]) % for win
end

% Compare the RMS values and find the best force model:
observation_set = ["night2_cr1: 0.6";
                   "night2_cr2: 1.2";
                   "night2_cr3: 1.8"];

best_rms_value = [night2_cr1.RMS.value;
                  night2_cr2.RMS.value;
                  night2_cr3.RMS.value];

cr_compare = table(observation_set,best_rms_value)

[cr_best_set,o] = min(cr_compare.best_rms_value);

fprintf('The best rms value we get from varying the solar reflectivity parameter is %0.8f\n', cr_best_set)
fprintf('which is from the %s observation set\n',cr_compare.observation_set(o))

format short
%%
save('study2_2_workspace.mat')

%% Save final project workspace

save('ENAE441_Project_workspace.mat')

%% Functions
function [obs1,obs2,obs3] = Get_observations(dataset, obs_indices)
%Get Observations
%   This function reads in dataset from one of the satellites and three desired 
%   observation indices and creates three structures for each of the three observations. 
%   Each observation structure includes the set of datetimes, lla, Right Ascension 
%   and Declinations. These three observations structures are used by the Gauss and Gibbs 
%   functions for initial orbit determination.

%   Written by: Ronak Chawla
obs1.time = dataset.datetime(obs_indices(1));
obs2.time = dataset.datetime(obs_indices(2));
obs3.time = dataset.datetime(obs_indices(3));

obs1.lla = latlonalt_deg(dataset.site_latitude_deg(obs_indices(1)), dataset.site_longitude_deg(obs_indices(1)), dataset.site_altitude_m(obs_indices(1)));
obs2.lla = latlonalt_deg(dataset.site_latitude_deg(obs_indices(2)), dataset.site_longitude_deg(obs_indices(2)), dataset.site_altitude_m(obs_indices(2)));
obs3.lla = latlonalt_deg(dataset.site_latitude_deg(obs_indices(3)), dataset.site_longitude_deg(obs_indices(3)), dataset.site_altitude_m(obs_indices(3)));

obs1.RA = dataset.right_ascension_deg(obs_indices(1));
obs2.RA = dataset.right_ascension_deg(obs_indices(2));
obs3.RA = dataset.right_ascension_deg(obs_indices(3));

obs1.Declination = dataset.declination_deg(obs_indices(1));
obs2.Declination = dataset.declination_deg(obs_indices(2));
obs3.Declination = dataset.declination_deg(obs_indices(3));

end


function [r1,r2,r3] = Gauss(ob1, ob2, ob3)
%Gauss Angles-Only function. 
%   Takes in three sets of angles (right ascension and declination) and 
%   their respective epochs to compute three position vectors r1, r2 and r3 
%   in units of km used by the Gibbs function for orbit determination.

%   Written by: Ronak Chawla and Prashanta Aryal

t1 = ob1.time; % datetime 1
t2 = ob2.time; % datetime 2
t3 = ob3.time; % datetime 3

% Right Ascension (alpha) and Declination of three observations in degrees
alpha1t = ob1.RA; 
delta1t = ob1.Declination; 

alpha2t = ob2.RA;  
delta2t = ob2.Declination;  

alpha3t = ob3.RA; 
delta3t = ob3.Declination; 

Tau1 = seconds(t1-t2);
Tau3 = seconds(t3-t2);

% Site Vectors from three observations in meters:

Rsv_1 = lla_to_eci(ob1.time, ob1.lla);
Rsv_2 = lla_to_eci(ob2.time, ob2.lla);
Rsv_3 = lla_to_eci(ob3.time, ob3.lla);

R = [Rsv_1 Rsv_2 Rsv_3]/1000; % in km

% Observation direction matrix

rho1_hat = [cosd(delta1t)*cosd(alpha1t); cosd(delta1t)*sind(alpha1t); sind(delta1t)];
rho2_hat = [cosd(delta2t)*cosd(alpha2t); cosd(delta2t)*sind(alpha2t); sind(delta2t)];
rho3_hat = [cosd(delta3t)*cosd(alpha3t); cosd(delta3t)*sind(alpha3t); sind(delta3t)];

L = [rho1_hat rho2_hat rho3_hat];

M = inv(L)*R;

a1 = Tau3/(Tau3 - Tau1);
a3 = -Tau1/(Tau3 - Tau1);

a1u = Tau3*((Tau3-Tau1)^2 - Tau3^2)/(6*(Tau3-Tau1));
a3u = -Tau1*((Tau3-Tau1)^2 - Tau1^2)/(6*(Tau3-Tau1));

A = M(2,1)*a1 - M(2,2) + M(2,3)*a3;
B = M(2,1)*a1u + M(2,3)*a3u;

R_2 = R(:,2);
E = dot(rho2_hat,R_2);

mu = 3.986*10^5;
coef3 = -mu^2*B^2;
coef2 = -2*mu*B*(A + E);
coef1 = -(A^2 + 2*A*E + norm(R_2)^2);

r2_roots = roots([1 0 coef1 0 0 coef2 0 0 coef3]);

 for i = 1:8
     
  if (r2_roots(i) > 6378) && (isreal(r2_roots(i)));
     r2 = r2_roots(i);
  end   
  
 end

u = mu/r2^3;

c1 = a1 + a1u*u;
c2 = -1;
c3 = a3 + a3u*u;

C = [-c1; -c2; -c3];

% rho1,2,3 in km:
rho = -(M*C)./C;

r1 = R(:,1) + rho(1)*rho1_hat;
r2 = R(:,2) + rho(2)*rho2_hat;
r3 = R(:,3) + rho(3)*rho3_hat;

end

function [v1, v2, v3] = Gibbs(r1, r2, r3)
%Gibbs method: 
%   This function takes in three position vectors in units of km of different
%   times and solves for the orbit by returning the velocity vectors in units 
%   of km/s corresponding to each position vector. The three points are coplanar.

%   Written by: Ronak Chawla

mu = 3.986*10^5;

D = cross(r2,r3) + cross(r3,r1) + cross(r1,r2);
N = norm(r1)*cross(r2,r3) + norm(r2)*cross(r3,r1) + norm(r3)*cross(r1,r2);
S = (norm(r1)-norm(r2))*r3 + (norm(r3)-norm(r1))*r2 + (norm(r2)-norm(r3))*r1;

B1 = cross(D,r1);
B2 = cross(D,r2);
B3 = cross(D,r3);

L = sqrt(mu/(norm(D)*norm(N)));

v1 = (L/norm(r1))*B1 + L*S;

v2 = (L/norm(r2))*B2 + L*S;

v3 = (L/norm(r3))*B3 + L*S;

end


function [RMS] = RMS(obs_dataset, pred_dataset)
%RMS Root Mean Square Function for angles.
%   [RMS] = RMS(obs_dataset, pred_dataset) returns the root mean square
%   value and residuals for azimuth and elevation.
%   
%   obs_dataset is a structure with fields azimuth_deg and elevation_deg
%   corresponding to an observed data set. pred_dataset is the exact same
%   except for a predicted dataset. Both variables should also have an
%   epoch field which is a column vector of datetimes. The length of both
%   obs_dataset and pred_dataset fields must all be the same length.
%
%   Written by: Ryan Quigley, Ronak Chawla, and Nico Lagendyk

azObs = obs_dataset.azimuth_deg;
elObs = obs_dataset.elevation_deg;

azPred = pred_dataset.azimuth_deg;
elPred = pred_dataset.elevation_deg;

timeObs = obs_dataset.epoch;

RMSsum = 0;

N = length(azObs);
azresidual = zeros(N,1);
elresidual = zeros(N,1);

    for i=1:N
        
        aztemp = mod(azObs(i) - azPred(i),360);
        azresidual(i) = min(360-aztemp,aztemp);
        eltemp = mod(elObs(i)-elPred(i),360);
        elresidual(i) = min(360-eltemp,eltemp);
        RMSsum = azresidual(i)^2 + elresidual(i)^2 + RMSsum;

    end

    RMS.value = sqrt((1/N)*RMSsum);
    RMS.azresidual = azresidual;
    RMS.elresidual = elresidual;
    
    plot(timeObs,azresidual)
    plot(timeObs,elresidual)
end

function [dataset] = ADD_aer(dataset, observer_lla)
%ADD_aer: Add Azimuth, Elevation and Range to structure
%
%   This function uses the position vectors at each time and an observer lla 
%   to convert the eci coordinates of a data into azimuth, elevation and range. 
%   The function uses a for loop to go through each time in the data set and
%   uses the Snag-App function eci_to_azelrn to obtain an azimuth, elevation 
%   and range for each time. Then, the function returns fields of azimuth
%   and elevation to the structure of the input dataset.

%   Written by: Ronak Chawla 

for i = 1:length(dataset.epoch)
    aer = eci_to_azelrn(dataset.epoch(i), dataset.position_m(i,:), observer_lla);
    dataset.azimuth_deg(i,1) = aer.azimuth_deg;
    dataset.elevation_deg(i,1) = aer.elevation_deg;
end

end

% Initial Orbit Determination function

function [output, night2_IOD_prop] = IOD(night1_dataset,night2_dataset,timespan,vector,name)
%IOD  Initial Orbit Determination
%
%   This function compiles all the steps necessary to go from a single set
%   of three indexes (vector) which specify rows in the night2_dataset.
%   These indexes are then used to pull observed values of the datetimes,
%   llas, right ascensions and declination using the Get_obervations
%   function.
%
%   Once the data is reformatted the Gauss function performs and angles
%   only calculation to find 3 position vectors. These vectors are then
%   plugged into the Gibbs funciton which solves for the corresponding
%   velocites to these positions thereby fully defining the orbit.
%
%   The three position velocity pairs are then converted to initial
%   estimates which are then propagated to times during the second night.
%   Then using the observed values from night2_dataset the RMS function is
%   used to compare the calculated. Then the rms values from the three
%   position velocity pairs are compared and the best (i.e. lowest
%   residual) is exported alongside all the individual rms values to the
%   output as fields in a structure. The RMS function also prints graphs
%   which are layered and exported as part of the output structure. 
%
%   The night2_IOD_prop are all the values from the propagation steps and
%   are also sent out with the associated values from night2_dataset which
%   allows for easy comparison.
%
%   Written by: Nico Lagendyk, Ronak Chawla, and Prashanta Aryal


    % Redefine variable names to make code easier to read
    ds1 = night1_dataset;
    ds2 = night2_dataset;
    ts = timespan;
    vec = vector;
    
    % Reformat observation data so it is easier to handle using
    % Get_obervations function
    [obs1,obs2,obs3] = Get_observations(ds1, vec);
    
    % Run Gauss code to determine position from 3 angles
    [r1, r2, r3] = Gauss(obs1,obs2,obs3); % r in units of km
    
    % Run Gibbs code to determine correlating velocities
    [v1, v2, v3] = Gibbs(r1, r2, r3); % v in units of km/s
    
    % Convert to m and m/s
    r1 = r1*1000; 
    r2 = r2*1000; 
    r3 = r3*1000;
    v1 = v1*1000; 
    v2 = v2*1000; 
    v3 = v3*1000;
    
    % Create 3 initial estimate pvts from calculated values
    initial_est1 = pvt(ds1.datetime(vec(1)),r1,v1);
    initial_est2 = pvt(ds1.datetime(vec(2)),r2,v2);
    initial_est3 = pvt(ds1.datetime(vec(3)),r3,v3);
    
    % Copy data from night 2 into separate structure
    night2_IOD_prop.obs.azimuth_deg = ds2.azimuth_deg(ts);
    night2_IOD_prop.obs.elevation_deg = ds2.elevation_deg(ts);
    night2_IOD_prop.obs.epoch = ds2.datetime(ts);
    
    force_model_4x4_1m2_cr1 = force_model_third_body(4, 4, 0, 0, 1, 1, 1000);
    observer_lla = latlonalt_deg(ds1.site_latitude_deg(1), ds1.site_longitude_deg(1), ds1.site_altitude_m(1));
    
    % Propagate inital estimate 1 to the same indicies of the second night
    night2_IOD_prop.pred1 = propagate_to_times(initial_est1, ds2.datetime(ts), force_model_4x4_1m2_cr1);
    
    for i = 1:length(night2_IOD_prop.pred1.epoch)
        aer1 = eci_to_azelrn(night2_IOD_prop.pred1.epoch(i), night2_IOD_prop.pred1.position_m(i,:), observer_lla);
        night2_IOD_prop.pred1.azimuth_deg(i,1) = aer1.azimuth_deg;
        night2_IOD_prop.pred1.elevation_deg(i,1) = aer1.elevation_deg;
    end
    
    % Propagate inital estimate 2 to the same indicies of the second night
    night2_IOD_prop.pred2 = propagate_to_times(initial_est2, ds2.datetime(ts), force_model_4x4_1m2_cr1);
    
    for i = 1:length(night2_IOD_prop.pred2.epoch)
        aer1 = eci_to_azelrn(night2_IOD_prop.pred2.epoch(i), night2_IOD_prop.pred2.position_m(i,:), observer_lla);
        night2_IOD_prop.pred2.azimuth_deg(i,1) = aer1.azimuth_deg;
        night2_IOD_prop.pred2.elevation_deg(i,1) = aer1.elevation_deg;
    end
    
    % Propagate inital estimate 3 to the same indicies of the second night
    night2_IOD_prop.pred3 = propagate_to_times(initial_est3, ds2.datetime(ts), force_model_4x4_1m2_cr1);
    
    for i = 1:length(night2_IOD_prop.pred3.epoch)
        aer1 = eci_to_azelrn(night2_IOD_prop.pred3.epoch(i), night2_IOD_prop.pred3.position_m(i,:), observer_lla);
        night2_IOD_prop.pred3.azimuth_deg(i,1) = aer1.azimuth_deg;
        night2_IOD_prop.pred3.elevation_deg(i,1) = aer1.elevation_deg;
    end
    
    % Setup figure for RMS plots
    output.rms_graph = figure;
    hold on
    
    % Determine RMS for each position velocity combination
    output.t1.rms = RMS(night2_IOD_prop.obs, night2_IOD_prop.pred1);
    output.t2.rms = RMS(night2_IOD_prop.obs, night2_IOD_prop.pred2);
    output.t3.rms = RMS(night2_IOD_prop.obs, night2_IOD_prop.pred3);
    
    title(sprintf('%s: Plot of Residuals over Time',name))
    xlabel('Date Time')
    ylabel('Residual (deg)')
    legend('t1 azimuth', 't1 elevation', 't2 azimuth', 't2 elevation', 't3 azimuth', 't3 elevation')
    
    if ismac
    saveas(gcf,[pwd sprintf('/ResidualGraphs_IOD/%s',name)]) % for mac
    elseif ispc
    saveas(gcf,[pwd sprintf('\\ResidualGraphs_IOD\\%s',name)]) % for win
    end
        
    
    % Return the initial estimates
    output.t1.initial_est = initial_est1;
    output.t2.initial_est = initial_est2;
    output.t3.initial_est = initial_est3;
    
    % Reformat vector sets into each sub-structure
    output.t1.posvel = [r1,v1];
    output.t2.posvel = [r2,v2];
    output.t3.posvel = [r3,v3];
    
    % Return Observation index into the structure
    output.t1.obs_index = vec(1);
    output.t2.obs_index = vec(2);
    output.t3.obs_index = vec(3);
    
    % Comparison of values
    rms_compare = [output.t1.rms.value, output.t2.rms.value, output.t3.rms.value];
    [output.best.rms, k] = min(rms_compare);
    output.best.index = vec(k);
    
end