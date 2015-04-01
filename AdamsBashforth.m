%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: [x, y] = AdamsBashforth(f_ode, xRange, yInitial, Tol, hMAx, hMin)
%
% Description:
% To approximate the solution of the initial-value problem
%       y = f (t, y), a ≤ t ≤ b, y(a) = α
% with local truncation error within a given tolerance. Using the
% Adams Bashforth 4th order Predictor and Adams-Moulton 3rd order 
% implicit Corrector.
% 
% Inputs: endpoints a, b; initial condition α; tolerance TOL; 
% maximum step size hmax; minimum step size hmin.
% Outputs: t,y vectors where at the ith step y(i,:) approximates y(ti)
%   
% Author:
% Pranav Mamidanna, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, y] = AdamsBashforth(f, tRange, yInitial, Tol, hMax, hMin)

    y(1,:) = yInitial;
    h = hMax;
    t(1,:) = tRange(1);
    FLAG = 1; LAST = 0;
    
    % Traditional Fourth Order Runge Kutta Method for Approximating
    % the first 3 values.
    [t, y] = RungeKutta4(f, 1, t, y, h);

    while(FLAG == 1)
        i = i+1;
        % Fourth Order Adams-Bashforth Predictor:
        yp(1,:) = y(i-1,:) + (h/24)*(55*(f(t(i-1),y(i-1))) - 59*(f(t(i-2),y(i-2)))...
                                   + 37*(f(t(i-3),y(i-3))) - 9*(f(t(i-4),y(i-4))));
        % Third Order Adams-Moulton Corrector:                      
        yc(1,:) = y(i-1,:) + (h/24)*(9*(f(t(i),yp)) + 19*(f(t(i-1),y(i-1)))...
                                   - 5*(f(t(i-2),y(i-2))) + 9*(f(t(i-3),y(i-3))));
        % Implementation of the Variable Step Size Algorithm
        sig = 19*abs(yp-yc)/(270*h);
        if sig <= Tol
            y(i,:) = yc(1,:);
            if LAST == 1
                FLAG = 0;
            else
                q = (Tol/(2*sig))^0.25;
                if q > 4
                    h = 4*h;
                else
                    h = h*q;
                end
                if h > hMax
                    h = hMax;
                end
                if t(i-1) + 4*h > xRange(2)
                    h = (xRange(2) - t(i-1))/4 ;
                    LAST = 1;
                end
                [t, y] = RungeKutta4(f, i-1, t, y, h);                
            end
        else
            q = (Tol/(2*sig))^0.25;
               if q < 0.1
                   h = 0.1*h;
               else
                   h = h*q;
               end
               if h < hMin
                   FLAG = 0;
                   print('hMin exceeded');
               else
                   [t, y] = RungeKutta4(f, i-1, t, y, h);
               end
        end
    end
    
end

function [t,y] = RungeKutta4(f, j, t, y, h) 
% 4th Order Runge Kutta Method to predict three consecutive
% values.
    for k = j:j+3
            t(j+1)  = t(j) + h;
            k1(1,:) = h*f(t(j), y(j,:));
            k2(1,:) = h*f(t(j) + h/2, y(j,:) + k1(1,:)/2);
            k3(1,:) = h*f(t(j) + h/2, y(j,:) + k2(1,:)/2);
            k4(1,:) = h*f(t(j+1), y(j,:) + k3(1,:));
            y(j+1,:) = y(j,:) + (k1(1,:) + 2*k2(1,:) + 2*k3(1,:) + k4(1,:)) * 1/6;
    end
end

% Example Usage:
% xa0 = 6; y0 = 0.029; zj0 = 0.0318;
% c0 = 2.0210e+13; c1 = 4.9687e+10;
% c2 = 15.6602; c4 = 133.7321; c3 = 20;
% 
% f = @(t,x) [xa0 - x(1) - c0*x(1)*(exp(-1/x(2)));
%             y0 - x(2) + c1*x(1)*(exp(-1/x(2))) - c2*(x(2) - x(3));
%             c3*(zj0 - x(3)) + c4*(x(2) - x(3))];
%         
% [t, y] = AdamsBashforth(f, [0 10], [1.36, 0.033, 0.0329], 10e-4, 0.01, 0.001);
% plot(t,y);        
