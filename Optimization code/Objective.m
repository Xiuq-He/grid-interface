function [objectiveQ,objectiveP] = Objective(alpha,Q,B,P,C)

%evaluation at current alpha iterate
t_ini_fcr = alpha(1);
t_a_fcr = alpha(2);
t_a_ffr = alpha(3);
t_d_ffr = alpha(4);
t_r_ffr = alpha(5);
x_ffr = alpha(6);
f_pod_low = alpha(7);
f_pod_high = alpha(8);
m_pod = alpha(9);
t_90_vc = alpha(10);
t_100_vc = alpha(11);


objectiveQ = trace(B(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)'*Q*B(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc));
objectiveP = trace(C(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)*P*C(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)');
end