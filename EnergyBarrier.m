function [Eb,Ebpp1,Ebpp2,Ebpplambda1,Ebpplambda2] = EnergyBarrier(x,l1R,l2R,l1D,l2D)
phi_1 = x(5);phi_2 = x(6);
eta_1 = phi_1;
phi_minus = min(phi_1,phi_2);
phi_plus = max(phi_1,phi_2);
interval_of_eta_2 = 100;

E_eta2 = zeros(interval_of_eta_2+2,2);
phi_34 = zeros(interval_of_eta_2+3,2);
phi_34(1,:) = phi_minus;
[s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4] = Para(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),l1R,l2R,l1D,l2D,0);

for k = 1:interval_of_eta_2+2
    eta_2 = phi_minus - (phi_plus - phi_minus)/(interval_of_eta_2 - 1) + (phi_plus - phi_minus)/(interval_of_eta_2 - 1)*(k - 1);
    myenergyfun = @(t)Energy_spr(eta_1,eta_2,t(1),t(2),s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4);
    options = optimset('MaxFunEvals',2e3,'MaxIter',2e3);
    [xx,pe] = fminsearch(@(t) myenergyfun(t),[0 0],options);
    phi_34(k+1,:) = xx;
    E_eta2(k,1) = eta_2;
    E_eta2(k,2) = pe;
end
%% find energy barrier
xx = E_eta2(:,1); yy = E_eta2(:,2);
kkk = find(xx < phi_plus & xx > phi_minus);
Eb = max(yy(kkk));

%% numerical stiffness (second derivative calculation via finite difference)
Ebpp1 = (E_eta2(3,2) + E_eta2(1,2) -2*E_eta2(2,2))/((phi_plus - phi_minus)/(interval_of_eta_2 - 1))^2;
Ebpp2 = (E_eta2(interval_of_eta_2,2) + E_eta2(interval_of_eta_2+2,2) -2*E_eta2(interval_of_eta_2+1,2))/((phi_plus - phi_minus)/(interval_of_eta_2 - 1))^2;

xxl = zeros(interval_of_eta_2+2,1);
xxtemp = E_eta2(:,1); 
for iiiii = 1:length(E_eta2)
    xxl(iiiii) = norm(Rot(phi_1)*s1 + Rot(xxtemp(iiiii))*s2 - Rot(phi_1)*t1 - Rot(xxtemp(iiiii))*t2);
end

Ebpplambda1 = (E_eta2(3,2) + E_eta2(1,2) - 2*E_eta2(2,2))/((xxl(3) - xxl(1))/2)^2;
Ebpplambda2 = (E_eta2(interval_of_eta_2,2) + E_eta2(interval_of_eta_2+2,2) - 2*E_eta2(interval_of_eta_2+1,2))/((xxl(interval_of_eta_2+2) - xxl(interval_of_eta_2))/2)^2;



end