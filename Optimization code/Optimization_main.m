clear all;
close all;
clc;

%Load input data
Optimization_data;

n_fcr = 4;
n_ffr = 6;
n_vq = 2;

%% projected gradient descent
steps= 10000; 

tic
[t_iterations,alpha_iterations,objQ_iterations,objP_iterations,maxQ_iterations,maxP_iterations]=Projected_GD(steps,alpha_0,dA_mf,dBB_mf,dCC_mf,grid_code_bounds,A_cl,B_cl,C_cl,BB_cl,CC_cl);
toc

%% evaluate optimization results
figure
plot(alpha_iterations(:,1),'LineWidth',2); hold on;
plot(alpha_iterations(:,2),'LineWidth',2);
plot(alpha_iterations(:,3),'LineWidth',2);
plot(alpha_iterations(:,4),'LineWidth',2);
plot(alpha_iterations(:,5),'LineWidth',2);
plot(alpha_iterations(:,6),'LineWidth',2);
plot(alpha_iterations(:,7),'LineWidth',2);
plot(alpha_iterations(:,8),'LineWidth',2);
plot(alpha_iterations(:,9),'LineWidth',2);
plot(alpha_iterations(:,10),'LineWidth',2);
plot(alpha_iterations(:,11),'LineWidth',2);
l = legend('$t_\mathrm{i}^\mathrm{fcr}$','$t_\mathrm{a}^\mathrm{fcr}$','$t_\mathrm{a}^\mathrm{ffr}$','$t_\mathrm{d}^\mathrm{ffr}$','$t_\mathrm{r}^\mathrm{ffr}$','$x_\mathrm{ffr}$','$f_\mathrm{low}^\mathrm{pod}$','$f_\mathrm{high}^\mathrm{pod}$','$m_\mathrm{pod}$','$t_\mathrm{90}^\mathrm{vq}$','$t_\mathrm{100}^\mathrm{vq}$');
set(l,'interpreter','latex','FontSize',12);


%% check feasibility of final alpha
final_feasibility = check_feasibility(alpha_iterations(end,:),grid_code_bounds);

if final_feasibility == 0
    fprintf('final alpha infeasible\n')    
end

if final_feasibility == 1
   fprintf('final alpha feasible\n')
end

%% evaluate alpha_iterations to get optimal T_des
s = tf('s');

t_ini_fcr = alpha_iterations(end,1);
t_a_fcr = alpha_iterations(end,2);
t_a_ffr = alpha_iterations(end,3);
t_d_ffr = alpha_iterations(end,4);
t_r_ffr = alpha_iterations(end,5);
x_ffr = alpha_iterations(end,6);
f_pod_low = alpha_iterations(end,7);
f_pod_high = alpha_iterations(end,8);
m_pod = alpha_iterations(end,9);
t_90_vc = alpha_iterations(end,10);
t_100_vc = alpha_iterations(end,11);


A_T_ss = A_T(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
B_T_ss = B_T(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
C_T_ss = C_T(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);

sys_T_ss = ss(A_T_ss,B_T_ss,C_T_ss,0);

if isstable(sys_T_ss) > 0
     fprintf('Tdes optimal ss stable\n')
end

if isstable(sys_T_ss) == 0
     fprintf('Tdes optimal ss unstable\n')
end


%% comparison with initial T_des

%FCR
t_ini_fcr_0 = 2;
t_a_fcr_0 = 30;

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
t_100_vc_0 = 60; 


A_T_ss_0 = A_T(t_ini_fcr_0,t_a_fcr_0,t_a_ffr_0,t_d_ffr_0,t_r_ffr_0,x_ffr_0,f_pod_low_0,f_pod_high_0,m_pod_0,t_90_vc_0,t_100_vc_0);
B_T_ss_0 = B_T(t_ini_fcr_0,t_a_fcr_0,t_a_ffr_0,t_d_ffr_0,t_r_ffr_0,x_ffr_0,f_pod_low_0,f_pod_high_0,m_pod_0,t_90_vc_0,t_100_vc_0);
C_T_ss_0 = C_T(t_ini_fcr_0,t_a_fcr_0,t_a_ffr_0,t_d_ffr_0,t_r_ffr_0,x_ffr_0,f_pod_low_0,f_pod_high_0,m_pod_0,t_90_vc_0,t_100_vc_0);

sys_T_ss_0 = ss(A_T_ss_0,B_T_ss_0,C_T_ss_0,0);

if isstable(sys_T_ss_0) > 0
     fprintf('Tdes initial ss stable\n')
end

if isstable(sys_T_ss_0) == 0
     fprintf('Tdes initial ss unstable\n')
end

figure
bode(sys_T_ss); hold on
bode(sys_T_ss_0);
legend('optimal','initial')

figure
step(-sys_T_ss,50); hold on
step(-sys_T_ss_0,50); 
legend('optimal','initial')

%% investigate closed-loop system response behavior
A_cl_ss = A_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
B_cl_ss = B_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
C_cl_ss = C_cl(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);

sys_cl_ss = ss(A_cl_ss,B_cl_ss,C_cl_ss,0);

 if isstable(sys_cl_ss) > 0
     fprintf('optimal closed-loop system stable\n')
 end

 if isstable(sys_cl_ss) == 0
     fprintf('optimal closed-loop system unstable\n')
 end
 

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

figure
step(-sys_cl_ss,50); hold on;
step(-sys_cl_ss0,50);
legend('optimal','initial')




