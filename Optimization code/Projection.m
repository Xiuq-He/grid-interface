function alpha_proj = Projection(alpha, gamma, gradient, grid_code_bounds)

fun = @(x) 0.5*norm((x - (alpha'-gamma.*gradient)),2)^2; %euclidean norm

x0 = alpha'; 

A = [-1 0 0 0 0 0 0 0 0 0 0;       % -t_i_fcr <= 0
     1 0 0 0 0 0 0 0 0 0 0;        %  t_i_fcr <= T_i_fcr_max
     1 -1 0 0 0 0 0 0 0 0 0;       %  t_i_fcr-t_a_fcr <= 0
     0 1 0 0 0 0 0 0 0 0 0;        %  t_a_fcr <= T_a_fcr_max
     1 -1 0 0 0 0 0 0 0 0 0;       %  t_i_fcr-t_a_fcr <= -|fcr_capacity|/Rp_ramp_max
     
     0 0 -1 0 0 0 0 0 0 0 0;       % -t_a_ffr <= 0
     0 0 1 0 0 0 0 0 0 0 0;        %  t_a_ffr <= T_a_ffr_max
     0 0 -1 0 0 0 0 0 0 0 0;       % -t_a_ffr <= -|ffr_capacity|/Rp_ramp_max
     0 0 -1 1 0 0 0 0 0 0 0;       %  t_d_ffr-t_a_ffr <= T_d_ffr_max
     0 0 1 -1 0 0 0 0 0 0 0;       %  t_a_ffr-t_d_ffr <= -T_d_ffr_min
     0 0 0 -1 1 0 0 0 0 0 0;       %  t_r_ffr-t_d_ffr <= T_r_ffr_max
     0 0 0 1 -1 0 0 0 0 0 0;       %  t_d_ffr-t_r_ffr <= -T_r_ffr_min
     0 0 0 0 0 -1 0 0 0 0 0;       % -x_ffr <= -1
     0 0 0 0 0 1 0 0 0 0 0;        %  x_ffr <= X_ffr_max
     0 0 0 0 0 1 0 0 0 0 0;        %  x_ffr <= Mp_peak_max/|ffr_capacity|
     
     0 0 0 0 0 0 -1 0 0 0 0;       % -f_pod_low <= -F_pod_min
     0 0 0 0 0 0 0 1 0 0 0;        %  f_pod_high <= F_pod_max
     0 0 0 0 0 0 1 -1 0 0 0;       %  f_pod_low - f_pod_high <= -9.2/T_a_ffr_max
     0 0 0 0 0 0 1 -1 0 0 0;       %  f_pod_low - f_pod_high <= -9.2/T_i_fcr_max
     
     0 0 0 0 0 0 0 0 0 -1 0;       % -t_90_vc <= 0
     0 0 0 0 0 0 0 0 0 1 0;        %  t_90_vc <= T_vc_max_90
     0 0 0 0 0 0 0 0 0 1 -1;       %  t_90_vc-t_100_vc <= 0
     0 0 0 0 0 0 0 0 0 0 1;        %  t_100_vc <= T_vc_max_100   
     0 0 0 0 0 0 0 0 0 -1 0;       % -t_90_vc <= -|q_capacity|/Rq_ramp_max 
     0 0 0 0 0 0 0 0 0 1 -1;       %  t_90_vc-t_100_vc <= -0.1*|q_capacity|/Rq_ramp_max
     ];

b = [-0.01;                                             % -t_i_fcr <= 0;
     grid_code_bounds(1);                               %  t_i_fcr <= T_i_fcr_max;
     0.01;                                              %  t_i_fcr-t_a_fcr <= 0
     grid_code_bounds(2);                               %  t_a_fcr <= T_a_fcr_max;
     grid_code_bounds(3)/grid_code_bounds(17);          %  t_i_fcr-t_a_fcr <= -|fcr_capacity|/Rp_ramp_max
     
     0.01;                                              % -t_a_ffr <= 0
     grid_code_bounds(4);                               %  t_a_ffr <= T_a_ffr_max
     grid_code_bounds(8)/grid_code_bounds(17);          % -t_a_ffr <= -|ffr_capacity|/Rp_ramp_max     
     grid_code_bounds(14);                              %  t_d_ffr-t_a_ffr <= T_d_ffr_max
     -grid_code_bounds(5);                              %  t_a_ffr-t_d_ffr <= -T_d_ffr_min
     grid_code_bounds(15);                              %  t_r_ffr-t_d_ffr <= T_r_ffr_max
     -grid_code_bounds(6);                              %  t_d_ffr-t_r_ffr <= -T_r_ffr_min
     -1.00001;                                          % -x_ffr <= -1
     grid_code_bounds(7);                               %  x_ffr <= X_ffr_max
     -grid_code_bounds(16)/grid_code_bounds(8);         %  x_ffr <= Mp_peak_max/|ffr_capacity|
     
     -grid_code_bounds(9);                              % -f_pod_low <= -F_pod_min
     grid_code_bounds(10);                              %  f_pod_high <= F_pod_max
     -9.2/grid_code_bounds(4);                          %  f_pod_low - f_pod_high <= -9.2/T_a_ffr_max
     -9.2/grid_code_bounds(1);                          %  f_pod_low - f_pod_high <= -9.2/T_i_fcr_max
     
     0.01;                                              % -t_90_vc <= 0
     grid_code_bounds(11);                              %  t_90_vc <= T_vc_max_90   
     0.01;                                              %  t_90_vc - t_100_vc <= 0
     grid_code_bounds(12);                              %  t_100_vc <= T_vc_max_100   
     grid_code_bounds(13)/grid_code_bounds(18);         % -t_90_vc <= -|q_capacity|/Rq_ramp_max 
     0.1*grid_code_bounds(13)/grid_code_bounds(18);     %  t_90_vc-t_100_vc <= -0.1*|q_capacity|/Rq_ramp_max
     ]; 



%nonlinear constraint for fmincon - with superposition constraints
    function [c,ceq] = nonlcon(x) 
        w_d = sqrt(4*((x(7)*x(8))/(x(8)-x(7))^2)-1);
        zeta = (x(8)-x(7))/2;
        c(1) = -grid_code_bounds(3)-grid_code_bounds(8)*x(6)+(2*x(9)/sqrt(w_d^2+1))*exp(-atan(w_d)/w_d) - grid_code_bounds(16);
        c(2) = -grid_code_bounds(3)-grid_code_bounds(8)*x(6)-((2*x(9)/sqrt(w_d^2+1))*exp(-atan(w_d)/w_d)) - grid_code_bounds(16);
        c(3) = -grid_code_bounds(3)/(x(2)-x(1))-grid_code_bounds(8)/x(3)+((2*x(9)/sqrt(w_d^2+1))*exp(-atan(w_d)/w_d))/(atan(w_d)/(zeta*w_d))-grid_code_bounds(17);
        c(4) = -grid_code_bounds(3)/(x(2)-x(1))-grid_code_bounds(8)/x(3)-((2*x(9)/sqrt(w_d^2+1))*exp(-atan(w_d)/w_d))/(atan(w_d)/(zeta*w_d))-grid_code_bounds(17);
        c(5) = 1-(4*x(7)*x(8))/(x(8)-x(7))^2;
        ceq = 0;
    end

%fmincon
alpha_proj_col = (fmincon(fun,x0,A,b,[],[],[],[],@nonlcon));
alpha_proj = alpha_proj_col';
 
end 