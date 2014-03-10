clear all
clc
syms alpha beta gamma real

% rotational matrices calculated in previous problem set
R_B1 = [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
R_12 = [cos(beta),0,sin(beta);0,1,0;-sin(beta),0,cos(beta)];
R_23 = [cos(gamma),0,sin(gamma);0,1,0;-sin(gamma),0,cos(gamma)];

% write down the 3x1 relative position vectors for link lengths l_i=1
r_3F_in3 = [0,0,-1]';
r_23_in2 = [0,0,-1]';
r_12_in1 = [0,0,-1]';
r_B1_inB = [0,1,0]';

% write down the homogeneous transformations
H_23 = [R_23,r_23_in2;0 0 0 1];
H_12 = [R_12,r_12_in1;0 0 0 1];
H_B1 = [R_B1,r_B1_inB;0 0 0 1];

% create the cumulative transformation matrix
H_B3 = H_B1*H_12*H_23; 

q = [alpha;beta;gamma];

% find the foot point position vector
r_BF_inB = H_B3(1:3,:)*[r_3F_in3;1];

% determine the foot point Jacobian J_BF_inB=d(r_BF_inB)/dq
J_BF_inB = [diff(r_BF_inB,q(1)) diff(r_BF_inB,q(2)) diff(r_BF_inB,q(3))];

% what generalized velocity dq do you have to apply in a configuration q = [0;60°;-120°]
% to lift the foot in vertical direction with v = [0;0;-1m/s];

qval = pi/180 * ([0;60;-120]);
Jval = eval(subs(J_BF_inB, q, qval));
dx = [0;0;-1];

dq = pinv(Jval) * dx;     

