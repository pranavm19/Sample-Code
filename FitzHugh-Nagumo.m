%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Course: BIO465.1x Neuronal Dynamics
% Instructor: Prof Wulfram Gurstner, EPFL
% 
% Description:
% Fitzhugh-Nagumo Model Simulation
% This 2-D model imitates spike generation by Hodgkin-Huxley type 
% "neurons". The equations are given by:
%       V' = V − (V^3)/3 − W + I
%       W' = 0.08(V + 0.7 − 0.8W)
%       Where,
%       V = Membrane Potential.
%       W = Recovery Variable (Capture Gating dynamics).
%       I = Injected Current.
% This code generates the Phase Plane Portrait of the system and 
% the Bifurcation Diagram as the injected stimulus current is 
% varied. 
%
% Author:
% Pranav Mamidanna, December 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Customary:
close all;
clear all;

% Global Variable: External Stimulus/Current Injected.
global i;
low = 0.2;
high = 1.2;
stepSize = 0.005;

for i = low:stepSize:high; 

f = @(w,v) [ v(1) - v(1).^3/3 - v(2) + i;
             0.7*(v(1) + 0.7 - 0.8*v(2)) ];
g = @(v) f(0,v);

% Finding the equilibria of the system at every "I" value.
% Deny fsolve ruining the command window.
options = optimset('Display','off');
fixedPoints = fsolve(g,[0 0], options);
Vs = fixedPoints(1);
Ws = fixedPoints(2);

% Jacobian of the linearised system.
% Eigenvalues and Coloring Scheme.
J = [ [1 - Vs^2, -1]; [0.7, -0.8*0.7]];
Lambda = eig(J);
if real(Lambda(1)) > 0; color = 'r.';
else color = 'b.';
end;

% Determining Nullclines. (dv/dt = 0, dw/dt = 0)
% Specify Range of Voltage.
v = -2:0.1:2;
wNull = (v + 0.7) / 0.8;
vNull = v - v.^3/3 + i;

% Phase Plane Portrait.
% Plot Nullclines.
subplot(1,2,1);
hold off;
plot(v, wNull, 'k-');
hold on;
plot(v, vNull, 'k-.');
% Plot trajectories generated from solving the ODEs.
[W, V] = ode45(f, [0 80], [Vs+0.2, Ws]);
plot(V(:,1), V(:,2), 'm-');
title(['Phase Portrait, I = ', num2str(i)]);
xlabel('V(t)'); ylabel('W(t)');
axis tight;

% Bifurcation Diagram.
subplot(1,2,2);
hold on;
[W, V] = ode45(f, [0 100], V(length(W),:) );
plot(i, min(V(:,1)), 'm.');
plot(i, max(V(:,1)), 'm.');
plot(i, Vs,color);
title('Bifurcation Diagram');
xlabel('I'); ylabel('V_e_q');
axis tight;
drawnow;

end;
