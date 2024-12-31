%% Simulation model of DC-coupled AC-output grid interface
% The model reproduces the result in the paper below.
% For more details, see the full text: https://arxiv.org/abs/2410.14912
% Developer: Xiuqiang He, ETH Zurich
% Last modification: 31 Dec. 2024
% MATLAB Version: R2022b

%% Simulation data for the IEEE 13 Node Test Feeder model
clear all;

% miles/km
mi2km = 1.60934;

% feet to km
ft2km = 0.0003048;

% microsiemens to Farads
ms2F = 1/2/pi/60*1e-6;


%% Configuration 601 - series reactance - ohm/mile

R_601 = [0.3465 0.1560 0.1580;0.1560 0.3375 0.1535;0.1580 0.1535 0.3414];
X_601 = [1.0179 0.5017 0.4236;0.5017 1.0478 0.3849;0.4236 0.3849 1.0348];

% charging susceptance - microsiemens/mile

B_601 = [6.2998 -1.9958 -1.2595;-1.9958 5.9597 -0.7417;-1.2595 -0.7417 5.6386];

% convert for SPS

R_601 = R_601/mi2km;

L_601 = X_601/mi2km/2/pi/60;

C_601 = B_601/mi2km*ms2F;


%% Configuration 602 - series reactance - ohm/mile

R_602 = [0.7526 0.1580 0.1560;0.1580 0.7475 0.1535;0.1560 0.1535 0.7436];
X_602 = [1.1814 0.4236 0.5017;0.4236 1.1983 0.3849;0.5017 0.3849 1.2112];

% charging susceptance - microsiemens/mile

B_602 = [5.6990 -1.0817 -1.6905;-1.0817 5.1795 -0.6588;-1.6905 -0.6588 5.4246];

% convert for SPS

R_602 = R_602/mi2km;

L_602 = X_602/mi2km/2/pi/60;

C_602 = B_602/mi2km*ms2F;


%% Configuration 603 - series reactance - ohm/mile

R_603 = [0 0 0;0 1.3294 0.2066;0 0.2066 1.3238];
X_603 = [0 0 0;0 1.3471 0.4591;0 0.4591 1.3569];

% charging susceptance - microsiemens/mile

B_603 = [0 0 0;0 4.7097 -0.8999;0 -0.8999 4.6658];

% convert for SPS

R_603 = R_603/mi2km;

L_603 = X_603/mi2km/2/pi/60;

C_603 = B_603/mi2km*ms2F;


%% Configuration 604 - series reactance - ohm/mile

R_604 = [1.3238 0 0.2066;0 0 0;0.2066 0 1.3294];
X_604 = [1.3569 0 0.4591;0 0 0;0.4591 0 1.3471];

% charging susceptance - microsiemens/mile

B_604 = [4.6658 0 -0.8999;0 0 0;-0.8999 0 4.7097];


% convert for SPS

R_604 = R_604/mi2km;

L_604 = X_604/mi2km/2/pi/60;

C_604 = B_604/mi2km*ms2F;


%% Configuration 605 - series reactance - ohm/mile

R_605 = [0 0 0;0 0 0;0 0 1.3292];
X_605 = [0 0 0;0 0 0;0 0 1.3475];

% charging susceptance - microsiemens/mile

B_605 = [0 0 0;0 0 0;0 0 4.5193];


% convert for SPS

R_605 = R_605/mi2km;

L_605 = X_605/mi2km/2/pi/60;

C_605 = B_605/mi2km*ms2F;

%% Configuration 606 - series reactance - ohm/mile

R_606 = [0.7982 0.3192 0.2849;0.3192 0.7891 0.3192;0.2849 0.3192 0.7982];
X_606 = [0.4463 0.0328 0.0143;0.0328 0.4041 0.0328;0.0143 0.0328 0.4463];
%X_606 = [0.4463 0.0328 -0.0143;0.0328 0.4041 0.0328;-0.0143 0.0328 0.4463];

% charging susceptance - microsiemens/mile

B_606 = [96.8897 -1e-6 -1e-6;-1e-6 96.8897 -1e-6;-1e-6 -1e-6 96.8897];

% convert for SPS

R_606 = R_606/mi2km;

L_606 = X_606/mi2km/2/pi/60;

C_606 = B_606/mi2km*ms2F;


%% Configuration 607 - series reactance - ohm/mile

R_607 = [1.3425 0 0;0 0 0;0 0 0];
X_607 = [0.5124 0 0;0 0 0;0 0 0];

% charging susceptance - microsiemens/mile

B_607 = [88.9912 0 0;0 0 0;0 0 0];

% convert for SPS

R_607 = R_607/mi2km;

L_607 = X_607/mi2km/2/pi/60;

C_607 = B_607/mi2km*ms2F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Below are the paremeters of the grid interface

% Frequency
f0 = 60;                        % Nominal frequency in Hz
w0 = 2*pi*f0;                   % Nominal angular frequency in rad/s
fw = 5e3;                       % PWM switching frequency in Hz

% Capacity
Stotal = 1e6;                   % Total capacity (VA) of AGI
SPV = Stotal * 0.7;             % PV capacity
SBess = Stotal * 0.3;           % BESS capacity
SLoad = Stotal * 0.1;           % Load
C = 80000e-6;                   % Size (uF) of the DC-link capacity
C_bc_PV = 5000e-6;              % Output capacitor of the PV boost converter
C_bc_ESS = 5000e-6;             % Output capacitor of the BESS boost converter
C_bc_SC = 5000e-6;              % Output capacitor of the SC boost converter
C_sum = C + C_bc_PV + C_bc_ESS + C_bc_SC;


% Voltage level
V1 = 4160;                      % Nominal AC voltage of the grid (rms ph-ph)
V2 = 400;                       % Nominal AC voltage of AGI (rms ph-ph)
Vdc_nom = 800;                  % DC voltage reference
Vac_nom = V2;
Vac_peak = Vac_nom*sqrt(2)/sqrt(3);   % AC voltage reference
Vdc_PV0 = 464.3;                % DC initial voltage of PV (deloaded, can
                                % be verified by clicking the "plot"
                                % command in the PV Array module.
Vdc_PV = 480;                   % DC nominal voltage of PV
Vdc_ESS = 480;                  % DC nominal voltage of ESS
Vdc_SCunit = 2.7;
Vdc_SC = Vdc_SCunit * 180;      % DC nominal voltage of SC
Vdc_SC0 = Vdc_SC;               % DC initial voltage of SC

% Per-unit base for high-voltage-side per-unit
Sbase = Stotal;
Vbase = V1;
Zbase = Vbase^2/Sbase;

% Per-unit base for DC AGI
Idcbase = Sbase/Vdc_nom;
Zdcbase = Vdc_nom^2/Sbase;
Cdcbase = 1/Zdcbase/w0;
C_pu = C_sum/Cdcbase;

% Base for supercap
C_SC = 20e3;                    % Size (F) of the supercapacitor set
C_SC_set = C_SC * Vdc_SCunit^2/Vdc_SC^2;
C_SC_eq = C_SC * Vdc_SCunit^2/(Vdc_nom^2 - (Vdc_nom*0.9)^2);
C_SC_pu = C_SC_eq/Cdcbase;
C_Inertia = 1/2 * C_SC_pu/w0;   % Inertia constant of the supercapacitor set

% Grid impedance
Rg = 0.07;
Zg = 0.2;

% Simulation setting
ts = 5e-6;                      % Simulation step size
tESSconnect = 0.1;              % Connection time of ESS
tSCconnect = 0.2;               % Connection time of SCap
tPVconnect = 0.3;               % Connection time of PV
tLoadconnect = 0.4;             % Connection time of PV
                                % Note: the start-up process is not the
                                % focus of this paper, so the start-up
                                % control is rough
tdist = 2;                      % Disturbance time
trcv = 3.5;
tend = 5;                       % End of the simulation
%% Specification and disaggregation of dynamic response
s = tf('s');
TJ = 5;                         % Desired inertial constant
D = 25;                         % Desired droop coefficient
T_des_vq = 0.2/(0.01 * s + 1);  % Desired voltage - var transfer function
                                % Please note the minus sign difference from
                                % Eq. (9) in the paper, similarly
                                % thereafter
T_des_vdcf = 10;
T_des_fvdc = 1/T_des_vdcf;
T_des_pf = (TJ*s + D);
T_des_fp = 1/T_des_pf;

T_SC = 0.2 * s;                 % Specification of SC response (0.0 - 0.4s)
T_SC1 = T_SC/(0.01*s + 1);                       % TF to be implemented (time constant 0.01 - 0.05 s)
T_agg = 1/(T_des_vdcf*T_des_fp) - 1/2*C_pu*s/w0 - T_SC; % - T_SC1 or - T_SC
T_PV = T_agg * 0.7/(0.5*s + 1);                   % Slow PV
T_ESS = T_agg - T_PV;                             % Fast ESS
T_ESS1 = T_ESS/(0.01*s + 1)/(0.01*s+1);           % TF to be implemented

% Re-aggregation for validation purpose
% (re-aggregated dynamics with the causalized desired control responses)
T_des = 1/(T_des_vdcf*T_des_fp);
T_des_agg = T_SC1 + T_PV + T_ESS1 + 1/2*C_pu*s/w0;
T_des_fp_agg_1 = 1/(T_des_vdcf*T_des_agg);

% Setpoint of the DC components (per-unit under Stotal)
Pset_PV = 0.5;
Pset_ESS = 0.0;
Pset_SC = 0.0;
Pset_load = -0.1; % Constant power load connected to the DC link

% Power limits
Plim_PV = 0.7;
Plim_ESS = 0.4;
Plim_SC = 0.3;

% Setpoints of the central inverter
Pset_inv = Pset_PV + Pset_ESS + Pset_load; % Only used for validation purpose
                                           % The AC-DC matching control
                                           % does not need a power
                                           % setpoint.
% Setpoints for GFM voltage magnitude control
Qset_inv = 0.0;
Vset_inv = 1.05;


%% Converter parameters

% Inverter parameters
Lf = 0.10;
Cf = 0.08;
Rf = 5e-4;
J = [0 -1; 1 0];

% Boost converter for PV
L_bc_PV = 4e-4;                 % L primary source side
R_bc_PV = 1e-3;                 % R primary source side
C_bc_PV = C_bc_PV;              % Output capacitor
R_line_PV = 1e-3;               % Line resistance

% Two-quarter converter for ESS
L_bc_ESS = 4e-4;                % L primary source side
R_bc_ESS = 1e-3;                % R primary source side
C_bc_ESS = C_bc_ESS;            % Output capacitor
R_line_ESS = 1e-3;              % Line resistance

% Two-quarter converter for SC
L_bc_SC = 4e-4;                 % L primary source side
R_bc_SC = 1e-3;                 % R primary source side
C_bc_SC = C_bc_SC;              % Output capacitor
R_line_SC = 1e-3;               % Line resistance

% PI gains for converter control
% Voltage control
Kp_v = 5;   % For PV
Ki_v = 500;
Kp_v1 = 50; % For SC and ESS
Ki_v1 = 1000;

% Current control
Kp_c = 10; 
Ki_c = 500;
Kp_c1 = 10;
Ki_c1 = 500;

% Filtering parameters
fc = 20;
fv = 100;
fc1 = 1000;

%% Current limit
Ilim = 1.2;