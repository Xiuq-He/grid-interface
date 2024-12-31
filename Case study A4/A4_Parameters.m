%% Simulation model of DC-coupled AC-output grid interface
% The model reproduces the result in the paper below.
% For more details, see the full text: https://arxiv.org/abs/2410.14912
% Developer: Xiuqiang He, ETH Zurich
% Last modification: 31 Dec. 2024
% MATLAB Version: R2022b

%% Parameters
clear all;

% Frequency
f0 = 60;                        % Nominal frequency in Hz
w0 = 2*pi*f0;                   % Nominal angular frequency in rad/s
fw = 5e3;                       % PWM switching frequency in Hz

% Capacity
Stotal = 1e6;                   % Total capacity (VA) of AGI
SPV = Stotal * 0.7;             % PV capacity
SPV1 = SPV * 0.7;
SPV2 = SPV * 0.3;
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

R1_line_PV = 1e-3;
R2_line_PV = 1e-3;
Vdc_nom1 = 480;
Vdc_nom2 = 480;
Vdc_PV1 = 200;
Vdc_PV2 = 200;
Vdc_PV0 = 480;
Vdc_PV10 = 232.1;
Vdc_PV20 = 232.1;

% Per-unit base for high-voltage-side per-unit
Sbase = Stotal;
Vbase = V1;
Zbase = Vbase^2/Sbase;

% Per-unit base for DC AGI
Idcbase = Sbase/Vdc_nom;
Idcbase1 = Sbase/Vdc_nom1;
Idcbase2 = Sbase/Vdc_nom2;
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
Rg = 0.04 * 2;                  % Grid resistance in per-unit (0.5 to 5 times)
Zg = 0.2 * 2;                   % Grid inductance in per-unit (0.5 to 5 times)

% Simulation setting
ts = 5e-6;                      % Simulation step size
tESSconnect = 0.1;              % Connection time of ESS
tSCconnect = 0.2;               % Connection time of SCap
tPVconnect = 0.3;               % Connection time of PV
tLoadconnect = 0.4;             % Connection time of PV
                                % Note: the start-up process is not the
                                % focus of this paper, so the start-up
                                % control is rough
tdist1 = 2;                     % Disturbance 1 (frequency dip) time
tdist2 = 4;                     % Disturbance 2 (voltage dip) time
tend = 5;                       % End time of the simulation


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
rho_des = 1;

T_SC = 0.2 * s;                 % Specification of SC response (0.0 - 0.4s)
T_SC1 = T_SC/(0.01*s + 1);                       % TF to be implemented (time constant 0.01 - 0.05 s)
T_agg = 1/(T_des_vdcf*T_des_fp) - 1/2*C_pu*s/w0 - T_SC; % - T_SC1 or - T_SC
T_PV = T_agg * 0.7/(0.5*s + 1);                   % Slow PV
T_PV1 = T_PV * 0.7;
T_PV2 = T_PV * 0.3;
T_ESS = T_agg - T_PV;                             % Fast ESS
T_ESS1 = T_ESS/(0.01*s + 1)/(0.01*s+1);           % TF to be implemented

T_PV0_vq = T_des_vq;
T_PV1_vq = T_des_vq;
T_PV2_vq = T_des_vq;
T_des_fvdc1 = T_des_fvdc;

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
Pset_PV1 = Pset_PV * 0.7;
Pset_PV2 = Pset_PV * 0.3;


% Power limits
Plim_PV = 0.7;
Plim_ESS = 0.4;
Plim_SC = 0.3;
Plim_PV1 = Plim_PV * 0.7;
Plim_PV2 = Plim_PV * 0.3;

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

% Inverter inner loop control parameters
fvc = 100; %200 50 oscillating
fvc_ff = 100; %100
kp_vc = 2; %2
ki_vc = 10; %10

fcc = 1000; %1000
fcc_ff = 50; %100
kp_cc = 0.8; %0.8
ki_cc = 20; %20