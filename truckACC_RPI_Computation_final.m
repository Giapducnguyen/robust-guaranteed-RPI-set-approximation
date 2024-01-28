% % New ACC Regulation Problem Formulation for Truck
close all; clear; clc;

%% Design Parameters

n_facet = 50;                           % number-of-facet threshold at which the approximation algorithm is activated
scaleRatio = [0.5; 0.5];               % hyperplane-location parameter (where to put the hyperplane)
nVR_maxIter = 1;                        % maximum number of volume reduction iterations
maxSegment = 7;                         % maximum segments of the line connecting a vertex and its projection on set Ek

% Best combination so far: 50 - 0.65 - 1 - 7

% Disturbance observer gain
l1 = 4.5;       %4.5;
L = [0, l1, 0];

%% Time Data
Ts = 0.1;                   % (s) sampling time

%% Model Data

% % first-order system representing input/output dynamics
% a_f = (K_L / (T_L*s + 1))*a_fdes
K_L = 1.0;                  % System gain
T_L = 0.45;                 % time constant

% % Constant time headway spacing policy
tau_h = 2.5;                % (s) nominal time headway
d0 = 5;                     % (m) stopping distance

% % State space model for car-following system
% X = [Delta d; Delta v; a_f]
% u = a_fdes
% d = v_p
% dX/dt = Phi*X + Pi*u + Gamma*d
Phi = [0, 1, -tau_h;
       0, 0, -1;
       0, 0, -1/T_L];
Pi =  [0; 0; K_L/T_L];
Gamma = [0; 1; 0];

[nx,nu] = size(Pi);
nd = size(Gamma, 2);

%% Modelling Discretization
nTaylor_iter = 100;        % number of summing iterations in Taylor expansion

% Car-following system
[Phid,Pid,Gammad] = getDiscreteModel(Phi,Pi,Gamma,Ts,nTaylor_iter);

%% System Constraints
Delta_d_min = -5;       % (m) minimum clearance error
Delta_d_max = 6;        % (m) maximum clearance error

Delta_v_min = -1;       % (m/s) minimum speed error
Delta_v_max = 0.9;      % (m/s) maximum speed error

a_min = -1.5;           % (m/s^2) minimum acceleration
a_max = 0.6;            % (m/s^2) maximum acceleration

%% RPI SET of Disturbance Estimation Error

% Acceleration increment set
delta_ap_min = -0.1;
delta_ap_max = 0.01;
set_delta_ap = Polyhedron('lb',delta_ap_min,'ub',delta_ap_max);
set_delta_d = set_delta_ap;

% Closed-loop disturbance estimation error matrix 
A_L = eye(nd) - L*Gammad;

% RPI set computation for disturbance estimation error
DO_RPI_outerApprox = true;

if ~DO_RPI_outerApprox
    figure
    hold on
    set_E1 = set_delta_d;
    set_E1.plot('color','r','alpha',0.02);
    n_max1 = 100;                            % Maximum number of iterations
    set_E1_compTime = nan(1, n_max1);
    set_E1_facetNumber = nan(1, n_max1);
    for n = 1:n_max1
        tic;
        FW1 = (A_L^n)*set_delta_d;
        set_E1 = set_E1 + FW1;
        set_E1.minHRep();
        set_E1.minVRep();
        set_E1_compTime(:,n) = toc;
        set_E1_facetNumber(:,n) = size(set_E1.getFacet(),1);
        set_E1.plot('color','r','alpha',0.02); drawnow;
    end
else
    
    epsilon1 = 5e-5;
    [s1,set_E1] = getRPIOuterApprox(epsilon1,A_L,set_delta_d,n_facet,scaleRatio,nVR_maxIter,maxSegment);
    
end
set_D = set_E1;

%% RPI SET of Actual-Nominal State Error

% car-following-system STATE set 
lb_X = [Delta_d_min; Delta_v_min; a_min];
ub_X = [Delta_d_max; Delta_v_max; a_max];
set_X = Polyhedron('lb',lb_X,'ub',ub_X);    

% car-following-system CONTROL set 
set_Uf = Polyhedron('lb',a_min,'ub',a_max); 

% car-following-system DISTURBANCE set
ap_min = -1.5;           % (m/s^2) minimum acceleration
ap_max = 0.6;            % (m/s^2) maximum acceleration
set_Up = Polyhedron('lb',ap_min,'ub',ap_max);
set_W = Gammad*set_D;

% Stabilization feedback gain
Q_fb = eye(nx);
R_fb = eye(nu);
[K_fb, ~, ~] = dlqr(Phid, Pid, Q_fb, R_fb, []);
Phid_cl = Phid - Pid*K_fb;  % closed-loop car-following-system state matrix
% spectral radius of the closed-loop car-following-system state matrix
spectral_radius = max(abs(eig(Phid_cl)));    

% RPI set Computation for car-following system
sys_RPI_outerApprox = true;

if ~sys_RPI_outerApprox
    figure
    hold on
    set_E2 = set_W;
    set_E2.plot('color','r','alpha',0.02);
    n_max2 = 20;                            % Maximum number of iterations
    set_E2_compTime = nan(1, n_max2);
    set_E2_facetNumber = nan(1, n_max2);
    set_E2_vertexNumber = nan(1, n_max2);
    for n = 1:n_max2
        fprintf('iter %d / %d\n',n,n_max2);
        tic;
        FW2 = (Phid_cl^n)*set_W;
        set_E2 = set_E2 + FW2;
        set_E2.minHRep();
        set_E2.minVRep();
        set_E2_compTime(:,n) = toc;
        set_E2_facetNumber(:,n) = size(set_E2.getFacet(),1);
        set_E2.plot('color','r','alpha',0.02); drawnow;
    end
else
    timeSetE = tic;
    epsilon2 = 1e-4;
    [s2,set_E2] = getRPIOuterApprox(epsilon2,Phid_cl,set_W,n_facet,scaleRatio,nVR_maxIter,maxSegment);
    comptime = toc(timeSetE);
end
set_X_tilde = set_E2;

%{
% % Save data
params_exact.L = L; params_exact.A_L = A_L;
params_exact.Ts = Ts;
params_exact.K_L = K_L;
params_exact.T_L = T_L;
params_exact.tau_h = tau_h;
params_exact.d0 = d0;
params_exact.Phi = Phi; params_exact.Pi = Pi; params_exact.Gamma = Gamma;
params_exact.nx = nx; params_exact.nu = nu; params_exact.nd = nd;
params_exact.A = A; params_exact.Bu = Bu;
params_exact.Phid = Phid; params_exact.Pid = Pid; params_exact.Gammad = Gammad;
params_exact.Ad = Ad; params_exact.Bud = Bud;
params_exact.Delta_d_min = Delta_d_min; params_exact.Delta_d_max = Delta_d_max;
params_exact.Delta_v_min = Delta_v_min; params_exact.Delta_v_max = Delta_v_max;
params_exact.a_min = a_min; params_exact.a_max = a_max;
params_exact.delta_ap_min = delta_ap_min; params_exact.delta_ap_max = delta_ap_max;
params_exact.ap_min = ap_min; params_exact.ap_max = ap_max; 
params_exact.set_D = set_D;
params_exact.set_X = set_X;
params_exact.set_Uf = set_Uf;
params_exact.set_W = set_W;
params_exact.set_X_tilde = set_X_tilde;
params_exact.Q_fb = Q_fb; params_exact.R_fb = R_fb; params_exact.K_fb = K_fb;
params_exact.Phid_cl = Phid_cl;
params_exact.n_facet = n_facet; params_exact.scaleRatio = scaleRatio;
save("params_exact.mat")
%}

% % ------------------|                    |-------------------------- % %
% % ------------------|  Function Helpers  |-------------------------- % %
% % ------------------|                    |-------------------------- % %

%% Discretization - ZOH
function [Ad,Bud,Bdd] = getDiscreteModel(Ac,Buc,Bdc,Ts,N_iter)
if det(Ac) ~= 0
    Ad = expm(Ac*ts);
    Bud = Ac\(Ad-eye(size(Ad)))*Buc;
    Bdd = Ac\(Ad-eye(size(Ad)))*Bdc;
else
    Ad = zeros(size(Ac));
    Bud = zeros(size(Buc));
    Bdd = zeros(size(Bdc));
    
    for k=0:N_iter
        Ad = Ad + (1/factorial(k))*(Ac*Ts)^k;
        Bud = Bud + (1/factorial(k+1))*Ac^(k)*Ts^(k+1)*Buc;
        Bdd = Bdd + (1/factorial(k+1))*Ac^(k)*Ts^(k+1)*Bdc;
    end
end
end

%% Function helper 1:
function [s,Fw_alpha] = getRPIOuterApprox(epsilon,Ak,set_W,n_facet,scaleRatio,nVR_iter,maxSegment)

[nx,~] = size(Ak);
Ms = 1000;
s = 0;
alpha = 1000;

mss = zeros(2*nx,1);

while(alpha > epsilon/(epsilon+Ms))
    s = s+1;
    alpha = max(set_W.support(Ak^s*set_W.A')./set_W.b);
    mss = mss + set_W.support([Ak^s, -Ak^s]);
    Ms = max(mss);
end
disp(s);
figure; hold on;
Fw = set_W; 
for i = 1:s-1
    % Fw.minHRep(); Fw.minVRep();
    
    set_temp = (Ak^i)*set_W;
    Fw = Fw + set_temp;
    Fw.minHRep(); Fw.minVRep();
    fprintf('iter %d / %d,  number of facets: %d \n',i,s, size(Fw.getFacet(),1));

    % If the number of facets of the current set > threshold
    % Proceed to approximate it

    if size(Fw.getFacet(),1) > n_facet
        % Ellipsoid-based polytope
        Fk_ebox = getEllipsoidPolytope(nVR_iter, scaleRatio, n_facet, Fw, maxSegment);
        
        % Extremum-based outer box
        Fk_obox = getExtremumOuterBox(Fw);
    
        % Intersection-based polytope 
        Fk_ibox = Fk_ebox & Fk_obox; 
        
        % update
        Fw = Fk_ibox;
        Fw.minHRep(); Fw.minVRep();
    end

    Fw.plot('color','red','alpha',0.01); drawnow limitrate nocallbacks
end
hold off;
Fw_alpha = ((1-alpha)^-1)*Fw;
end