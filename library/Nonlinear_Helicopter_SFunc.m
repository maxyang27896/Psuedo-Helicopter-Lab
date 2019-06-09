function [sys,x0,str,ts] = Nonlinear_Helicopter_SFunc(t,x,u,flag,varargin)
% S-Function style nonlinear model for the lab helicopter setup
% Equations based on Quanser document 'Dynamic Equations for the 3-DOF Helicopter'
% Version 1 (21/11/2016)


persistent DERX 

switch flag,
  % Initialization 
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes(varargin);

  % Derivatives 
  case 1,
   sys = DERX;
   
  % Update 
  case 2,
    sys=[];

  % Outputs 
  case 3,
   [sys,DERX] = NonlinearHelicopter(t,x,u,varargin);
   

  % GetTimeOfNextVarHit 
  case 4,
    sys=[];

  % Terminate 
  case 9,
    sys=[];

  % Unexpected flags 
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end


function [sys,x0,str,ts]=mdlInitializeSizes(varargin)

sizes = simsizes;

sizes.NumContStates  = 6;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 2;
sizes.DirFeedthrough =  1;
sizes.NumSampleTimes =  1;   % at least one sample time is needed

sys = simsizes(sizes);

% initialize the initial conditions

% user supplied IC
if nargin > 0,
  x0  = varargin{1}{2};
else % default IC
  x0  = [0;0;0;0;0;0];
end;

% str is always an empty matrix
str = [];

% initialize the array of sample times
ts  = [0 0];

function [Y,DERX] = NonlinearHelicopter(t,X,U,varargin)

% user supplied parameters
if nargin > 3,
    p = varargin{1}{1};
    % Propeller force-thrust constant found experimentally (N/V)
    Kf = p(1);
    % Mass of the helicopter body (kg)
    mh = p(2);
    % Mass of counter-weight (kg)
    mw = p(3);
    % Mass of front propeller assembly = motor + shield + propeller + body (kg)
    mf = p(4);
    % Mass of back propeller assembly = motor + shield + propeller + body (kg)
    mb = p(5);
    % Distance between pitch pivot and each motor (m)
    Lh = p(6);
    % Distance between elevation pivot to helicopter body (m)
    La = p(7);
    % Distance between elevation pivot to counter-weight (m)
    Lw = p(8);
    % Gravitational Constant (m/s^2)
    g =  p(9);    
else %default parameters 
    % Propeller force-thrust constant found experimentally (N/V)
    Kf = 0.1188;
    % Mass of the helicopter body (kg)
    mh = 1.308;
    % Mass of counter-weight (kg)
    mw = 1.924;
    % Mass of front propeller assembly = motor + shield + propeller + body (kg)
    mf = mh / 2;
    % Mass of back propeller assembly = motor + shield + propeller + body (kg)
    mb = mh / 2;
    % Distance between pitch pivot and each motor (m)
    Lh = 7.0 * 0.0254;
    % Distance between elevation pivot to helicopter body (m)
    La = 26.0 * 0.0254;
    % Distance between elevation pivot to counter-weight (m)
    Lw = 18.5 * 0.0254;
    % Gravitational Constant (m/s^2)
    g = 9.81;    
end;

%Nonlinear System
epsilon_dot=X(4);
rho_dot=X(5);
lambda_dot=X(6);
epsilon_ddot = -cos(X(2)) ^ 2 * Lh * (mf - mb) * La * sin(X(2)) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) * X(4) ^ 2 + (-0.2e1 * Lh * cos(X(2)) * (cos(X(1)) * La * (mf - mb) * cos(X(2)) ^ 2 - cos(X(1)) * La * (mf - mb) + sin(X(2)) * sin(X(1)) * Lh * (mf + mb)) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) * X(6) - 0.2e1 * (mf + mb) * cos(X(2)) * Lh ^ 2 * sin(X(2)) * X(5) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw)) * X(4) + ((sin(X(2)) * La * (mf - mb) * cos(X(1)) - 0.2e1 * sin(X(1)) * Lh * (mf + mb)) * Lh * cos(X(1)) * cos(X(2)) ^ 2 - 0.2e1 * sin(X(2)) * La * Lh * (mf - mb) * cos(X(1)) ^ 2 - sin(X(1)) * ((mf + mb) * La ^ 2 - Lh ^ 2 * mb - Lh ^ 2 * mf + Lw ^ 2 * mw) * cos(X(1)) + sin(X(2)) * La * Lh * (mf - mb)) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) * X(6) ^ 2 - 0.2e1 * Lh * (cos(X(1)) * Lh * (mf + mb) * cos(X(2)) ^ 2 - cos(X(1)) * Lh * (mf + mb) - sin(X(2)) * sin(X(1)) * La * (mf - mb)) * X(5) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) * X(6) + Lh * (mf - mb) * La * sin(X(2)) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) * X(5) ^ 2 + (((g * mb + g * mf + Kf * (U(1) + U(2))) * La - Lw * mw * g) * ((mf - mb) ^ 2 * La ^ 2 + Lh ^ 2 * (mf + mb) ^ 2) * cos(X(2)) ^ 3 + (g * (mw * Lw * (mf - mb) ^ 2 * La ^ 2 + ((-0.4e1 * Lh ^ 2 * mf + Lw ^ 2 * mw) * mb ^ 2 + (-0.4e1 * Lh ^ 2 * mf ^ 2 - 0.2e1 * Lw ^ 2 * mw * mf) * mb + mf ^ 2 * mw * Lw ^ 2) * La + mw * Lh ^ 2 * Lw * (mf + mb) ^ 2) * cos(X(1)) + ((g * mb + g * mf + Kf * (U(1) + U(2))) * La - Lw * mw * g) * ((mf - mb) ^ 2 * La ^ 2 + Lh ^ 2 * (mf + mb) ^ 2)) * cos(X(1)) * cos(X(2)) ^ 2 + (Kf * (-U(1) + U(2)) * (mf - mb) * ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) * La * cos(X(1)) - (sin(X(1)) * Kf * Lh * (-U(1) + U(2)) * sin(X(2)) + (g * mb + g * mf + Kf * (U(1) + U(2))) * La - Lw * mw * g) * ((mf - mb) ^ 2 * La ^ 2 + Lh ^ 2 * (mf + mb) ^ 2)) * cos(X(2)) + (0.4e1 * La ^ 2 * mb * mf + mw * Lw ^ 2 * (mf + mb)) * cos(X(1)) * (-g * ((mf + mb) * La - Lw * mw) * cos(X(1)) + g * sin(X(1)) * Lh * (mf - mb) * sin(X(2)) + (g * mb + g * mf + Kf * (U(1) + U(2))) * La - Lw * mw * g)) / (0.4e1 * La ^ 2 * mb * mf + mw * Lw ^ 2 * (mf + mb)) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) / cos(X(1));
rho_ddot = cos(X(2)) * (((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) * sin(X(2)) * cos(X(1)) + cos(X(2)) ^ 2 * sin(X(1)) * La * Lh * (mf - mb)) / cos(X(1)) / ((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) * X(4) ^ 2 + (-0.2e1 * (-cos(X(2)) ^ 2 * ((mf + mb) * La ^ 2 + Lw ^ 2 * mw) * cos(X(1)) ^ 2 + cos(X(2)) ^ 2 * sin(X(2)) * sin(X(1)) * La * Lh * (mf - mb) * cos(X(1)) - Lh ^ 2 * (mf + mb) * cos(X(2)) ^ 2 + (mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) / cos(X(1)) / ((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) * X(6) + 0.2e1 * (mf + mb) * cos(X(2)) ^ 2 * Lh ^ 2 * sin(X(1)) * X(5) / cos(X(1)) / ((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw)) * X(4) - cos(X(2)) * (((-mb - mf) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) * sin(X(2)) * cos(X(1)) ^ 3 + sin(X(1)) * La * Lh * (cos(X(2)) ^ 2 - 0.2e1) * (mf - mb) * cos(X(1)) ^ 2 + 0.2e1 * sin(X(2)) * Lh ^ 2 * (mf + mb) * cos(X(1)) + sin(X(1)) * La * Lh * (mf - mb)) / cos(X(1)) / ((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) * X(6) ^ 2 - 0.2e1 * Lh * cos(X(2)) * (-La * (mf - mb) * cos(X(1)) ^ 2 + sin(X(2)) * sin(X(1)) * Lh * (mf + mb) * cos(X(1)) + La * (mf - mb)) * X(5) / cos(X(1)) / ((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) * X(6) - cos(X(2)) * Lh * sin(X(1)) * (mf - mb) * La / cos(X(1)) / ((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) * X(5) ^ 2 + (-g * cos(X(2)) * (mf - mb) * ((-0.4e1 * La * mb * mf + mw * Lw * (mf + mb)) * La * Lh ^ 2 + Lw * ((mf + mb) * La ^ 2 + Lw ^ 2 * mw) * (La + Lw) * mw) * cos(X(1)) ^ 3 + (-Lh ^ 2 * Kf * (-U(1) + U(2)) * ((mf - mb) ^ 2 * La ^ 2 + Lh ^ 2 * (mf + mb) ^ 2) * cos(X(2)) ^ 2 + (g * Lh * sin(X(1)) * ((mf + mb) * (-0.4e1 * La * mb * mf + mw * Lw * (mf + mb)) * Lh ^ 2 + Lw * La * mw * (mf - mb) ^ 2 * (La + Lw)) * sin(X(2)) - ((g * mb + g * mf + Kf * (U(1) + U(2))) * La - Lw * mw * g) * (mf - mb) * ((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) * La) * cos(X(2)) - Kf * (-U(1) + U(2)) * ((-mb - mf) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) * ((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw)) * cos(X(1)) ^ 2 + (-((g * mb + g * mf + Kf * (U(1) + U(2))) * La - Lw * mw * g) * (mf - mb) * ((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) * La * cos(X(2)) ^ 2 + Lh * (((g * mb + g * mf + Kf * (U(1) + U(2))) * La - Lw * mw * g) * sin(X(1)) * ((mf - mb) ^ 2 * La ^ 2 + Lh ^ 2 * (mf + mb) ^ 2) * sin(X(2)) - g * Lh * (0.4e1 * La ^ 2 * mb * mf + mw * Lw ^ 2 * (mf + mb)) * (mf - mb)) * cos(X(2)) + (0.2e1 * sin(X(1)) * Kf * Lh * (-U(1) + U(2)) * sin(X(2)) + (g * mb + g * mf + Kf * (U(1) + U(2))) * La - Lw * mw * g) * (mf - mb) * ((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) * La) * cos(X(1)) + Lh * (((g * mb + g * mf + Kf * (U(1) + U(2))) * La - Lw * mw * g) * sin(X(1)) * sin(X(2)) + Kf * Lh * (-U(1) + U(2))) * (((mf - mb) ^ 2 * La ^ 2 + Lh ^ 2 * (mf + mb) ^ 2) * cos(X(2)) ^ 2 - (mf + mb) * ((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw))) / Lh / (0.4e1 * La ^ 2 * mb * mf + mw * Lw ^ 2 * (mf + mb)) / ((mf + mb) * Lh ^ 2 + (mf + mb) * La ^ 2 + Lw ^ 2 * mw) / cos(X(1)) ^ 2;
lambda_ddot = -cos(X(2)) ^ 3 * Lh * (mf - mb) * La / cos(X(1)) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) * X(4) ^ 2 + (0.2e1 * (Lh * (sin(X(2)) * La * (mf - mb) * cos(X(1)) - sin(X(1)) * Lh * (mf + mb)) * cos(X(2)) ^ 2 + sin(X(1)) * ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw)) / cos(X(1)) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) * X(6) - 0.2e1 * (mf + mb) * cos(X(2)) ^ 2 * Lh ^ 2 * X(5) / cos(X(1)) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw)) * X(4) + cos(X(2)) * Lh * (La * (mf - mb) * cos(X(1)) ^ 2 * cos(X(2)) ^ 2 + 0.2e1 * sin(X(2)) * sin(X(1)) * Lh * (mf + mb) * cos(X(1)) - 0.2e1 * (mf - mb) * (cos(X(1)) ^ 2 - 0.1e1 / 0.2e1) * La) / cos(X(1)) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) * X(6) ^ 2 + 0.2e1 * (cos(X(1)) * Lh * (mf + mb) * sin(X(2)) + sin(X(1)) * La * (mf - mb)) * cos(X(2)) * Lh * X(5) / cos(X(1)) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) * X(6) + cos(X(2)) * Lh * (mf - mb) * La / cos(X(1)) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) * X(5) ^ 2 + (-(((g * mb + g * mf + Kf * (U(1) + U(2))) * La - Lw * mw * g) * sin(X(2)) + sin(X(1)) * Kf * Lh * (-U(1) + U(2))) * ((mf - mb) ^ 2 * La ^ 2 + Lh ^ 2 * (mf + mb) ^ 2) * cos(X(2)) ^ 2 - ((g * (mw * Lw * (mf - mb) ^ 2 * La ^ 2 + ((-0.4e1 * Lh ^ 2 * mf + Lw ^ 2 * mw) * mb ^ 2 + (-0.4e1 * Lh ^ 2 * mf ^ 2 - 0.2e1 * Lw ^ 2 * mw * mf) * mb + mf ^ 2 * mw * Lw ^ 2) * La + mw * Lh ^ 2 * Lw * (mf + mb) ^ 2) * cos(X(1)) + ((g * mb + g * mf + Kf * (U(1) + U(2))) * La - Lw * mw * g) * ((mf - mb) ^ 2 * La ^ 2 + Lh ^ 2 * (mf + mb) ^ 2)) * sin(X(2)) - g * Lh * sin(X(1)) * (0.4e1 * La ^ 2 * mb * mf + mw * Lw ^ 2 * (mf + mb)) * (mf - mb)) * cos(X(1)) * cos(X(2)) + ((-Kf * La * (-U(1) + U(2)) * (mf - mb) * cos(X(1)) + (mf + mb) * ((g * mb + g * mf + Kf * (U(1) + U(2))) * La - Lw * mw * g)) * sin(X(2)) + sin(X(1)) * Kf * Lh * (-U(1) + U(2)) * (mf + mb)) * ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw)) / (0.4e1 * La ^ 2 * mb * mf + mw * Lw ^ 2 * (mf + mb)) / ((mf + mb) * La ^ 2 + Lh ^ 2 * mb + Lh ^ 2 * mf + Lw ^ 2 * mw) / cos(X(1)) ^ 2;

%State Derivative
DERX=[epsilon_dot;rho_dot;lambda_dot;epsilon_ddot;rho_ddot;lambda_ddot];

%Output
Y=[X(1);X(2);X(3);epsilon_dot;rho_dot;lambda_dot];
