function gamma = Gamma_computation(gradient,k)


stepsize = 0.3*[0.02;     %t_ini_fcr 
            0.02;         %t_a_fcr
            0.4;          %t_a_ffr 
            0.1;          %t_d_ffr
            0.1;          %t_r_ffr
            0.1;          %x_ffr
            0.05;         %f_pod_low
            0.01;         %f_pod_high
            0.05;         %m_pod
            0.1;          %t_vq_90 
            0.1];         %t_vq_100 


%additional conditions   
if k > 5000
    gamma = abs(stepsize./gradient)*(5000/k)^10;
else
    gamma = abs(stepsize./gradient);
end




end
