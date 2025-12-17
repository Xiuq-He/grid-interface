function [t,alpha,objQ,objP,maxQeq,maxPeq] = Projected_GD(steps, alpha_0, dA, dBB, dCC, grid_code_bounds, A, B, C, BB, CC)

alpha(1,:)= alpha_0;

for k = 1:steps
    [gradient,Q,P,maxQ,maxP] = Gradient_computation(alpha(k,:),dA,dBB,dCC,A,BB,CC);
    [objQ(k),objP(k)] = Objective(alpha(k,:),Q,B,P,C);
    maxQeq(k) = maxQ;
    maxPeq(k) = maxP;
    gamma = Gamma_computation(gradient,k);
    gradient_gamma = gamma.*gradient
    alpha(k+1,:) = Projection(alpha(k,:),gamma,gradient,grid_code_bounds);
    
    %termination condition 
    if k > 100
    if (norm(alpha(k+1,:)-alpha(k-99,:),2) <= 1*10^(-3)) && (norm(alpha(k+1,:)-alpha(k-9,:),2) <= 1*10^(-3))
    fprintf('alpha converged\n')
    return
    end
    end
    
    t(k) = k-1;
end

end