function [gradient,Q,P,maxQ,maxP] = Gradient_computation(alpha, dA_mf, dBB_mf, dCC_mf, A, BB, CC)

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


%solve Lyapunov equations
P = lyap(A(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc),BB(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc));
Q = lyap(A(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)',CC(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc));

%check P and Q matrices
maxP=max(max(abs(A(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)*P+P*A(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)'+BB(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc))));
maxQ=max(max(abs(A(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)'*Q+Q*A(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc)+CC(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc))));

%gradient computation
dA1 = dA_mf{1}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dA2 = dA_mf{2}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dA3 = dA_mf{3}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dA4 = dA_mf{4}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dA5 = dA_mf{5}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dA6 = dA_mf{6}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dA7 = dA_mf{7}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dA8 = dA_mf{8}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dA9 = dA_mf{9}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dA10 = dA_mf{10}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dA11 = dA_mf{11}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);



dBB1 = dBB_mf{1}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dBB2 = dBB_mf{2}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dBB3 = dBB_mf{3}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dBB4 = dBB_mf{4}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dBB5 = dBB_mf{5}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dBB6 = dBB_mf{6}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dBB7 = dBB_mf{7}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dBB8 = dBB_mf{8}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dBB9 = dBB_mf{9}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dBB10 = dBB_mf{10}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dBB11 = dBB_mf{11}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);


dCC1 = dCC_mf{1}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dCC2 = dCC_mf{2}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dCC3 = dCC_mf{3}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dCC4 = dCC_mf{4}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dCC5 = dCC_mf{5}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dCC6 = dCC_mf{6}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dCC7 = dCC_mf{7}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dCC8 = dCC_mf{8}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dCC9 = dCC_mf{9}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dCC10 = dCC_mf{10}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);
dCC11 = dCC_mf{11}(t_ini_fcr,t_a_fcr,t_a_ffr,t_d_ffr,t_r_ffr,x_ffr,f_pod_low,f_pod_high,m_pod,t_90_vc,t_100_vc);


gradient(1,:)= 2*trace(dA1*P*Q)+trace(dBB1*Q)+trace(P*dCC1);
gradient(2,:)= 2*trace(dA2*P*Q)+trace(dBB2*Q)+trace(P*dCC2);
gradient(3,:)= 2*trace(dA3*P*Q)+trace(dBB3*Q)+trace(P*dCC3);
gradient(4,:)= 2*trace(dA4*P*Q)+trace(dBB4*Q)+trace(P*dCC4);
gradient(5,:)= 2*trace(dA5*P*Q)+trace(dBB5*Q)+trace(P*dCC5);
gradient(6,:)= 2*trace(dA6*P*Q)+trace(dBB6*Q)+trace(P*dCC6);
gradient(7,:)= 2*trace(dA7*P*Q)+trace(dBB7*Q)+trace(P*dCC7);
gradient(8,:)= 2*trace(dA8*P*Q)+trace(dBB8*Q)+trace(P*dCC8);
gradient(9,:)= 2*trace(dA9*P*Q)+trace(dBB9*Q)+trace(P*dCC9);
gradient(10,:)= 2*trace(dA10*P*Q)+trace(dBB10*Q)+trace(P*dCC10);
gradient(11,:)= 2*trace(dA11*P*Q)+trace(dBB11*Q)+trace(P*dCC11);




end
