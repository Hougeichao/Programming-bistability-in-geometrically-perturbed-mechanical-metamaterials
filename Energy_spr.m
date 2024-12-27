function [PE] = Energy_spr(eta_1,eta_2,eta_3,eta_4,s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4)

R1 = Rot(eta_1); R2 = Rot(eta_2); R3 = Rot(eta_3); R4 = Rot(eta_4);
PE = 2*norm(-1/2*(R1*s1 + R2*s2 + R3*s3 + R4*s4) + 3/8*(R1*t1 + R2*t2 + R3*t3 + R4*t4) + 1/8*(R1*u1 + R2*u2 + R3*u3 + R4*u4) - 1/4*(R1*v1 + R2*v2 + R3*v3 + R4*v4)).^2 + ...
     2*norm(1/4*(R1*s1 + R2*s2 + R3*s3 + R4*s4) - 1/8*(R1*t1 + R2*t2 + R3*t3 + R4*t4) - 3/8*(R1*u1 + R2*u2 + R3*u3 + R4*u4) + 1/2*(R1*v1 + R2*v2 + R3*v3 + R4*v4)).^2 + ...
     2*norm(-1/8*(R1*s1 + R2*s2 + R3*s3 + R4*s4) + 1/2*(R1*t1 + R2*t2 + R3*t3 + R4*t4) + 1/4*(R1*u1 + R2*u2 + R3*u3 + R4*u4) - 3/8*(R1*v1 + R2*v2 + R3*v3 + R4*v4)).^2 + ...
     2*norm(3/8*(R1*s1 + R2*s2 + R3*s3 + R4*s4) - 1/4*(R1*t1 + R2*t2 + R3*t3 + R4*t4) - 1/2*(R1*u1 + R2*u2 + R3*u3 + R4*u4) + 1/8*(R1*v1 + R2*v2 + R3*v3 + R4*v4)).^2;

end

