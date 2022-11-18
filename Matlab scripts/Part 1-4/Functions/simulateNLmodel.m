function newstate = simulateNLmodel(x0,u,Ts)
% Simulates the non-linear model without that super slow, stupid simulink
% model.
tspan = [0 Ts];
[t,x] = ode45(@(t,x) vdPol_continuous(t,x,u), tspan, x0);
newstate = x(end,:)';

% continuous vd pol equation with fixed u.
    function dxdt = vdPol_continuous(t,x,u)
        x1 = x(1); x2 = x(2);
        mu = 2;
        dxdt = [x2;
            mu*(1-x1^2)*x2-x1+u];
    end

end