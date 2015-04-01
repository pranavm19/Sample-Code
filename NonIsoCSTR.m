%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Deterministic Chaos in a Non-Isothermal CSTR.
% 
% Description:
% A Stirred Tank Reactor for Propylene Glycol Production
% is shown to exhibit chaotic response. The differential
% equations (in their dimensionless form) describing the 
% dynamics are:
%       Xa' = Xa0 - Xa - c0.Xa.exp(-1/Y)
%       Xb' = Xb0 - Xa - c01.Xa.exp(-1/Y)
%       Xc' = -Xc - c02.Xa.exp(-1/Y)
%       Xm' = Xm0 - Xm
%       Y'  = Y0  - Y  + c1.Xa.exp(-1/Y) - c2.(Y-Z)
%       Z'  = c3.(Zj0-Z) + c4.(Y-Z)
%       where,
%       Xa,Xb,Xc,Xd - Dimensionless Concentrations
%       Y, Z - Dimensionless Reactor and Jacket Temp. resp.  
% 
% The system is reduced to a 3-dimensional system containing 
% the variables Xa, Y, Z. This system is explored for existence 
% of oscillating regimes, unstability etc.
%
% This particular code produces graphs, depicting:
% 1. Chaotic Response
% 2. Eigenvalue Analysis
% 3. Chaos Map, which is the parameter values for sinusoidal
%    variation of jacket flow rate, for which chaos is observed.
%
% Author:
% Pranav Mamidanna, March 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializing System Parameters.
y0 = 0.0285; zj0 = 0.0318;
c0 = 2.0210e+13; c1 = 4.9687e+10;
c2 = 15.6602; c4 = 133.7321;
% Sinusoidal Perturbation Parameters.
w = 1.9; c3b = 22.76; cm3 = 9;

% Simulating two very close initial conditions:
for xa0 = [5, 5.0001];
    % Function Handle:
    f = @(t,x) [xa0 - x(1) - c0*x(1)*(exp(-1/x(2)));
            y0 - x(2) + c1*x(1)*(exp(-1/x(2))) - c2*(x(2) - x(3));
            (c3b + cm3*sin(w*t))*(zj0 - x(3)) + c4*(x(2) - x(3))];
    if xa0 == 5 ; [T1, X1] = ode45(f, [0, 100], [1.36, 0.0337, 0.0335]);
    else [T2, X2] = ode45(f, [0, 100], [1.36, 0.0337, 0.0335]);
    end;
end

%%%%
% Chaotic response of the system.
%%%%
subplot(2,2,1)
axis tight
title('System Response')
plot(T1, X1(:,2)*100 - X2(:,2)*100)
subplot(2,2,2)
title('Phase Plane Portrait')
axis tight
plot(X1(:,1), X1(:,2))

%%%%
% Eigenvalue Analysis for the System
%%%%
Lambda = zeros(3,(8-5)/0.1);
i = 1;
xa0 = 5:0.1:8;
for xa0 = [5:0.1:8];
    % Function Handle:
    f = @(t,x) [xa0 - x(1) - c0*x(1)*(exp(-1/x(2)));
            y0 - x(2) + c1*x(1)*(exp(-1/x(2))) - c2*(x(2) - x(3));
            c3b*(zj0 - x(3)) + c4*(x(2) - x(3))];
    g = @(x) f(0,x);
    options = optimset('Display','off');
    fixedPoints = fsolve(g, [0 0 0], options);
    xs = fixedPoints(1); ys = fixedPoints(2);
    Jacobian = [[ (-1 -c0*(exp(-1/ys))), (c0*(xs*((exp(-1/ys))))), 0];
            [(c1*(exp(-1/ys))), (-1-c2+c1*(xs*((exp(-1/ys))))), c2];
            [0, c4, (-c3b-c4)]];
    Lambda(:,i) = eig(Jacobian); i=i+1;
end
subplot(2,2,3)
axis tight
plot(Lambda(1,:)); hold on
axis tight
plot(Lambda(2,:),'k'); hold off

%%%%
% Plotting the Chaotic Map.
% cm3, w are the parameter for sinusoidal variation of flow in cooling jacket.
%%%%
A = []; B = [];
for cm3 = 4:0.5:10;                     % Ideal Stepsize = 0.01 or less
    for w = 1:0.2:3                     % Then a detailed map is obtained.
        clear X1 T1 X2 T2;
        for xa0 = [5, 5.0001]
            % Function Handle:
            f = @(t,x) [xa0 - x(1) - c0*x(1)*(exp(-1/x(2)));
                        y0 - x(2) + c1*x(1)*(exp(-1/x(2))) - c2*(x(2) - x(3));
                        (c3b + cm3*sin(w*t))*(zj0 - x(3)) + c4*(x(2) - x(3))];
            if xa0 == 5 ; [T1, X1] = ode45(f, [0:0.01:100], [1.36, 0.0337, 0.0335]);
            else [T2, X2] = ode45(f, [0:0.01:100], [1.36, 0.0337, 0.0335]);
            end;
        end
        Y = abs(X2(:,2)*100-X1(:,2)*100);
        if max(Y) > 0.001
            A = [A cm3];
            B = [B w];
        end
    end
end
subplot(2,2,4)
axis tight
plot(B,A,'r.')
