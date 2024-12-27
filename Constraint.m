function [g,h] = Constraint(x,l1R,l2R,l1D,l2D)
fac = 1; 
epsilon = 0.01*fac*fac;
epsilon_2 = 0.01;
theta_1 = x(5); theta_2 = x(6); theta_3 = x(7); theta_4 = x(8);
R1 = Rot(theta_1);R2 = Rot(theta_2);R3 = Rot(theta_3);R4 = Rot(theta_4);
W = [0,-1;1,0];
[s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4] = Para(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),l1R,l2R,l1D,l2D,0);

%% reference slit 
g(1) = -dot(W*s1,s2) + epsilon;
g(2) = -dot(W*s2,s3) + epsilon;
g(3) = -dot(W*s3,s4) + epsilon;
g(4) = -dot(W*s4,s1) + epsilon;
g(5) = -dot(t1,W*t2) + epsilon;
g(6) = -dot(t2,W*t3) + epsilon;
g(7) = -dot(t3,W*t4) + epsilon;
g(8) = -dot(t4,W*t1) + epsilon;
g(9) = -dot(u1,W*u2) + epsilon;
g(10) = -dot(u2,W*u3) + epsilon;
g(11) = -dot(u3,W*u4) + epsilon;
g(12) = -dot(u4,W*u1) + epsilon;
g(13) = -dot(W*v1,v2) + epsilon;
g(14) = -dot(W*v2,v3) + epsilon;
g(15) = -dot(W*v3,v4) + epsilon;
g(16) = -dot(W*v4,v1) + epsilon;

%% panels
g(17) = -dot(W*t1,s1) + epsilon;
g(18) = -dot(W*v1,t1) + epsilon;
g(19) = -dot(W*u1,v1) + epsilon;
g(20) = -dot(W*s1,u1) + epsilon;
g(21) = -dot(W*u2,s2) + epsilon;
g(22) = -dot(W*v2,u2) + epsilon;
g(23) = -dot(W*t2,v2) + epsilon;
g(24) = -dot(W*s2,t2) + epsilon;
g(25) = -dot(W*t3,s3) + epsilon;
g(26) = -dot(W*v3,t3) + epsilon;
g(27) = -dot(W*u3,v3) + epsilon;
g(28) = -dot(W*s3,u3) + epsilon;
g(29) = -dot(W*u4,s4) + epsilon;
g(30) = -dot(W*v4,u4) + epsilon;
g(31) = -dot(W*t4,v4) + epsilon;
g(32) = -dot(W*s4,t4) + epsilon;

%% deformed slits
g(33) = -dot(W*R1*s1,R2*s2) + epsilon;
g(34) = -dot(W*R2*s2,R3*s3) + epsilon;
g(35) = -dot(W*R3*s3,R4*s4) + epsilon;
g(36) = -dot(W*R4*s4,R1*s1) + epsilon;
g(37) = -dot(R1*t1,W*R2*t2) + epsilon;
g(38) = -dot(R2*t2,W*R3*t3) + epsilon;
g(39) = -dot(R3*t3,W*R4*t4) + epsilon;
g(40) = -dot(R4*t4,W*R1*t1) + epsilon;
g(41) = -dot(R1*u1,W*R2*u2) + epsilon;
g(42) = -dot(R2*u2,W*R3*u3) + epsilon;
g(43) = -dot(R3*u3,W*R4*u4) + epsilon;
g(44) = -dot(R4*u4,W*R1*u1) + epsilon;
g(45) = -dot(W*R1*v1,R2*v2) + epsilon;
g(46) = -dot(W*R2*v2,R3*v3) + epsilon;
g(47) = -dot(W*R3*v3,R4*v4) + epsilon;
g(48) = -dot(W*R4*v4,R1*v1) + epsilon;

g(49) = -((det(R1 - R2))^2 + (det(R2 - R3))^2 + (det(R3 - R4))^2 + (det(R4 - R1))^2) + epsilon_2;
h = [];

end