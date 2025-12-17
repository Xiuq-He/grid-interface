%% define symbolic parameters 
syms t_ini_fcr t_a_fcr real
syms t_a_ffr t_d_ffr t_r_ffr x_ffr real
syms f_pod_low f_pod_high m_pod real
syms t_90_vc t_100_vc real

%assumptions
assume(t_ini_fcr,'positive'); hold on;
assume(t_a_fcr,'positive');
assume(t_a_ffr,'positive'); 
assume(t_d_ffr,'positive');
assume(t_r_ffr,'positive');
assume(x_ffr,'positive');
assume(f_pod_low,'positive');
assume(f_pod_high,'positive');
assume(t_90_vc,'positive');
assume(t_100_vc,'positive');
close all;

%% initial parameters
%initial parameter vector alpha_0 (needs to be feasible)
%FCR
t_ini_fcr_0 = 0.01;
t_a_fcr_0 = 20;

%FFR
t_a_ffr_0 = 2;
t_d_ffr_0 = t_a_ffr_0+8;
t_r_ffr_0 = t_d_ffr_0+10;
x_ffr_0 = 1.00001;

%POD
f_pod_low_0 = 2*pi*1;
f_pod_high_0 = 2*pi*2;
m_pod_0 = -5; 

%VQ
t_90_vc_0 = 5;
t_100_vc_0 = 20; 

alpha_0 = [t_ini_fcr_0,t_a_fcr_0,t_a_ffr_0,t_d_ffr_0,t_r_ffr_0,x_ffr_0,f_pod_low_0,f_pod_high_0,m_pod_0,t_90_vc_0,t_100_vc_0];


%% parameter bounds given by the grid code requirements and device-level limitations
%FCR
T_i_fcr_max = 2;
T_a_fcr_max = 30;
fcr_capacity = -1/0.05;

%FFR
T_a_ffr_max = 2;
T_d_ffr_min = 8;
T_r_ffr_min = 10;
X_ffr_max = 1.35;
ffr_capacity = -1/0.04;

%POD                                            
F_pod_min = 2*pi*0.1;                  
F_pod_max = 2*pi*3;         

%VQ
T_vc_max_90 = 5;
T_vc_max_100 = 60;
q_capacity = -1/0.04; 

%device-level limits
T_d_ffr_max = 20;
T_r_ffr_max = 20;
Mp_peak_max = 70;
Rp_ramp_max = 111.3636;
Rq_ramp_max = 150;

grid_code_bounds = [T_i_fcr_max,
    T_a_fcr_max,
    fcr_capacity,
    T_a_ffr_max,
    T_d_ffr_min,
    T_r_ffr_min,
    X_ffr_max,
    ffr_capacity,
    F_pod_min,
    F_pod_max,
    T_vc_max_90,
    T_vc_max_100,
    q_capacity,
    T_d_ffr_max,
    T_r_ffr_max,
    Mp_peak_max,
    Rp_ramp_max,
    Rq_ramp_max
    ];

initial_feasibility = check_feasibility(alpha_0,grid_code_bounds);

if initial_feasibility == 0
    fprintf('initial alpha infeasible\n')
    return;
    
end

if initial_feasibility == 1
   fprintf('initial alpha feasible\n')
end


%% build closed-loop system matrices
%ss matrices of grid equivalence estimation A_G, B_G, Bh_G, C_G, where
% [\Delta \dot f_pll
%  \Delta f_pll      
%  \Delta \int f_pll  =  Z_g * [\Delta p
%  \Delta v                     \Delta q]
%  \Delta \int v]                               
load('grid_equivalence_nominal.mat')

%ss matrices of Tdes A_T, B_T, C_T, where
% [\Delta p        = Tdes *  [\Delta f_pll
%  \Delta q]                  \Delta v]
load('parametric_Tdes_mf.mat')

%specify performance output z = [\Delta\dot f_pll,\Delta f_pll,\Delta\int f_pll,\Delta v,\Delta\int v]
m_dot_f_pll = 0;
m_f_pll = 2*10^(1);
m_int_f_pll = 2*10^(1);  
m_v = 0;
m_int_v = 7*10^(-1);
M_weight = diag([m_dot_f_pll,m_f_pll,m_int_f_pll,m_v,m_int_v]);
M_half = M_weight^(1/2);



%parametric closed-loop system matrices
A_cl = @(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc) [A_G, B_G*C_T(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc); B_T(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)*[C_G(2,:);C_G(4,:)], A_T(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)];
B_cl = @(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc) [Bh_G; zeros(size(A_T(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc),1),2)];
C_cl = @(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc) [M_half*C_G, zeros(5,size(A_T(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc),1))];



%% check initial closed-loop system stability
A_cl_ss0 = A_cl(t_ini_fcr_0,t_a_fcr_0,t_a_ffr_0,t_d_ffr_0,t_r_ffr_0,x_ffr_0,f_pod_low_0,f_pod_high_0,m_pod_0,t_90_vc_0,t_100_vc_0);
B_cl_ss0 = B_cl(t_ini_fcr_0,t_a_fcr_0,t_a_ffr_0,t_d_ffr_0,t_r_ffr_0,x_ffr_0,f_pod_low_0,f_pod_high_0,m_pod_0,t_90_vc_0,t_100_vc_0);
C_cl_ss0 = C_cl(t_ini_fcr_0,t_a_fcr_0,t_a_ffr_0,t_d_ffr_0,t_r_ffr_0,x_ffr_0,f_pod_low_0,f_pod_high_0,m_pod_0,t_90_vc_0,t_100_vc_0);

sys_cl_ss0 = (ss(A_cl_ss0,B_cl_ss0,C_cl_ss0,0));

 if isstable(sys_cl_ss0) > 0
     fprintf('initial closed-loop system stable\n')
 end

 if isstable(sys_cl_ss0) == 0
     fprintf('initial closed-loop system unstable\n')
 end
 
%prescaling
[sys_cl_ss0_ps,info_ps] = prescale(sys_cl_ss0);
T_l = diag(info_ps.SL);
T_r = diag(info_ps.SR);

%re-scaled system matrices
A_cl = @(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)T_l*A_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)*T_r;
B_cl = @(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)T_l*B_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
C_cl = @(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)C_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)*T_r;

BB_cl = @(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc) B_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)*B_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)';
CC_cl = @(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc) C_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)'*C_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);



%% symbolic precomputations
%symbolic parameter vector alpha
alpha_parameters = [t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc];

%symbolic gradient computation
n_parameters = length(alpha_parameters);

for k = 1:n_parameters
    dA_sym = diff(A_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc),alpha_parameters(k));
    dA_mf{k} = matlabFunction(dA_sym,'Vars',[t_ini_fcr t_a_fcr t_a_ffr t_d_ffr t_r_ffr x_ffr f_pod_low f_pod_high m_pod t_90_vc t_100_vc]);

    dBB_sym = diff(BB_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc),alpha_parameters(k));
    dBB_mf{k} = matlabFunction(dBB_sym,'Vars',[t_ini_fcr t_a_fcr t_a_ffr t_d_ffr t_r_ffr x_ffr f_pod_low f_pod_high m_pod t_90_vc t_100_vc]);

    dCC_sym = diff(CC_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc),alpha_parameters(k));
    dCC_mf{k} = matlabFunction(dCC_sym,'Vars',[t_ini_fcr t_a_fcr t_a_ffr t_d_ffr t_r_ffr x_ffr f_pod_low f_pod_high m_pod t_90_vc t_100_vc]);
end 

