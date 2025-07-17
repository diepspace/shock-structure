%% 
close all
clear
clc


%% Computational parameters (user defined)
% Guess itermax and dx depending on the assumed shock thickness
itermax=15000;  % Maximum marching steps
dx=-1e-9;       % Marching direction/distance, m
epsi=1e-8;      % Deviation at starting condition to begin march


%% Flow parameters (user defined)
% Calorically perfect gas
gamma=1.4;
R=287;                % m^2/s^2-K
cp=R*gamma/(gamma-1);
mu0=1.716e-5;         % N-s/m^2
n_mu=0.666;
k0=0.0241;            % W/m-K
n_k=0.81;
T0=273;               % K


%% Shock upstream parameters (station 1) (user defined)
M1=1.1;
T1=290;   % K
p1=101e3; % Pa


%% Shock upstream conditions (station 1) (computed)
rho1=p1/(R*T1);
a1=sqrt(gamma*R*T1);
u1=M1*a1;
mu1=mu0*(T1/T0)^n_mu;
h1=cp*T1;
h01=h1+u1.^2/2;


%% Conserved quantities across the shock (here computed using upstream values)
A=rho1*u1;
B=rho1*u1^2+p1;
C=rho1*u1*h01;


%% Shock downstream conditions (station 2)
% Computed from conserved quantities A, B, C
u2=(B-sqrt(B^2-4*A*(1-R/2/cp)*(R*C/cp)))/(2*A*(1-R/2/cp));
rho2=A/u2;
T2=1/cp*(C/A-u2^2/2);
p2=rho2*R*T2;


%% Midpoint velocity used for coordinate shift calculation
u_mid=(u1+u2)/2;


%% Euler method (marching scheme) to evaluate equations derived from NS equations across shock
x=0;
u=(1+epsi)*u2; % Marching scheme requires slight gradient to begin
T=(1-epsi)*T2;

for iter=1:itermax
    % Recasted Navier-Stokes equations (dudx and dTdx)
    p_old=A/u(end)*R*T(end);
    mu_old=mu0*(T(end)/T0)^n_mu;
    k_old=k0*(T(end)/T0)^n_k;
    h_old=cp*T(end);

    dudx_old=3/4/mu_old*(A*u(end)+p_old-B); % Momentum equation
    u_new=u(end)+dudx_old*dx; % Euler method

    dTdx_old=1/k_old*(A*(h_old-u(end)^2/2)+(B-p_old)*u(end)-C); % Energy equation
    T_new=T(end)+dTdx_old*dx;

    x_new=x(end)+dx; % March in dx direction

    % Add latest step to entire vector (x, u, T span entire shock width)
    x=[x;x_new];
    u=[u;u_new];
    T=[T;T_new];

    % Compute value to shift all x-values such that x-x0=0 corresponts to middle of shock thickness
    ind_mid=find(u>=u_mid);
    if ~isempty(ind_mid)
        ind_mid=ind_mid(1);
        x_0=x(ind_mid); % Value to shift x
    else
        x_0=0;
    end

    % Integration stopping criteria
    if abs((u(end)-u1)/u1)<=1e-3&&abs(x(end))>=abs(2*x_0)
        break
    end
end


%% Remaining flow properties after integration
rho=A./u;
p=rho*R.*T;
ds=cp*log(T/T1)-R*log(p/p1);


%% Shock Thickness (Poggie Textbook Method)
ind_xa=find((u1-u)/(u1-u2)<=0.01);
    if ~isempty(ind_xa)
        ind_xa=ind_xa(1);
    end
ind_xb=find((u1-u)/(u1-u2)<=0.99);
    if ~isempty(ind_xb)
        ind_xb=ind_xb(1);
    end
xt=abs(x(ind_xa)-x(ind_xb));

Re_xt=rho1*u1*xt/mu1; % Nondim thickness in terms of Reynolds number
lambda=mu1/(0.67*rho1*a1); % Mean free path
xt_lambda=xt/lambda; % Nondim thickness in terms of mean free path


%% Plot internal shock structure
x_nondim=rho1*u1*(x-x_0)/mu1; % Nondim distance along shock width

figure
hold on
plot(x_nondim,u/u1,'DisplayName','u/u_1')
plot(x_nondim,T/T1,'DisplayName','T/T_1')
plot(x_nondim,p/p1,'DisplayName','p/p_1')
plot(x_nondim,rho/rho1,'DisplayName','\rho/\rho_1')
plot(x_nondim,ds/R,'DisplayName','\Deltas/R')
hold off
grid on
legend(gca,'show')
xlabel('\rho_{1}u_{1}(x-x_0)/\mu_1')

title_string=sprintf('Internal Shock Structure (M_1 = %1.1f)',M1);
title(title_string)