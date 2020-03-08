function [xp]=solver_ode45(t,x,red_A,red_B,b0,hyd_force,vane_pass_force,phi,phi_dot,phi_double_dot)
% xp = red_A*x + real(red_B*b0*(phi_dot^2-phi_double_dot)*exp(1i*phi) + red_B*hyd_force*exp(1i*phi) + red_B*vane_pass_force*exp(1i*4*phi));
xp = red_A*x + real(red_B*b0*(phi_dot^2-phi_double_dot)*exp(1i*phi));
end
