clear; clc;

q_0 = 1e-05; % m^3/s
V = 1.75e-3; % m^3
beta_e = 1e+009; % Pa
p_3 = 1e5; % Pa
D_m = 6.67E-05/(2*pi); % m^3/rad
T_L = 300; % Nm
J_t = 0.1; % kg m^2
B_m = 0.3; % Nm/rad

s = tf('s');

A = beta_e/ (V*s);
B = D_m / (J_t*s + B_m);

C = A / ( 1 + A*B*D_m );
D = ( A*B*T_L - p_3 ) / A;

omega_h = sqrt( (beta_e*D_m^2) / (V*J_t) )
sigma_h = (B_m*V*omega_h) / (2*beta_e*D_m^2)
bode(C)