function f = Obj(x,l1R,l2R,l1D,l2D)

[Eb] = EnergyBarrier(x,l1R,l2R,l1D,l2D); 
f = (Eb - 0.003)^2;

% [s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4,Module_1_xi,Module_1_lambda] = Para(s11,s12,s21,s22,phi_1,phi_2,phi_3,phi_4,l1R,l2R,l1D,l2D,0);
% C1 = 500;C2 = 1;
% f = (C1*(Eb - 0.001)^2 + C2*(Module_1_lambda - 0.46)^2);
end