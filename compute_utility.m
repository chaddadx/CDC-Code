function [U] = compute_utility(gamma1, gamma2,T, k, lambda, alpha12, alpha21)

 % Set initial conditions
    X1 = zeros(T,1);
    X2 = zeros(T,1);
   
    X1(1) = 1;
    X2(1) = 1;
    
    for t = 2:T
        X1(t) = k * ( gamma1 * X1(t-1) + alpha21 * X2(t-1) );
        X2(t) = k * ( gamma2 * X2(t-1) + alpha12 * X1(t-1) );
    end
    
    cumX1 = sum(X1);
    cumX2 = sum(X2);
    
   
    U1 = cumX1 - lambda * gamma1;
    U2 =  cumX2 - lambda * gamma2;
   
    % Combine utilities into a vector
    U = [U1, U2];
end







