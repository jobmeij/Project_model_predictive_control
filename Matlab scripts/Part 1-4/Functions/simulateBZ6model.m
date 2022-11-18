function newstate = simulateBZ6model(x0,u,Ts)
% Simulates the non-linear model without that super slow, stupid simulink
% model.
tspan = [0 Ts];
[t,x] = ode45(@(t,x) BZ6_continuous(t,x,u), tspan, x0);
newstate = x(end,:)';

% continuous vd pol equation with fixed u.
    function dxdt = BZ6_continuous(t,x,u)
        x1 = x(1); x2 = x(2);
        
        dxdt = [x2;
            -0.33*exp(-x1)*x1-1.1*x2+u];
    end

end