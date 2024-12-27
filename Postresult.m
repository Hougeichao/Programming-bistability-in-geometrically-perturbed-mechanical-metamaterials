%% Run Postresult.m after implementing optimization to plot energy results %%
tic
phi_1 = x(5);phi_2 = x(6);
eta_1 = phi_1;
phi_minus = min(phi_1,phi_2);
phi_plus = max(phi_1,phi_2);
interval_of_eta_2 = 5000;

E_eta2 = zeros(interval_of_eta_2,2);
[s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4] = Para(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),l1R,l2R,l1D,l2D,0);
phi_34 = zeros(interval_of_eta_2,2);

for k = 1:interval_of_eta_2
    eta_2 = phi_minus - 0.2 + (phi_plus + 0.2 - phi_minus + 0.2)/(interval_of_eta_2 - 1)*(k - 1);
    myenergyfun = @(t)Energy_spr(eta_1,eta_2,t(1),t(2),s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4);
    options = optimset('MaxFunEvals',2e3,'MaxIter',2e3);
    [xxxx,pe] = fminsearch(@(t) myenergyfun(t),[0 0],options);
    phi_34(k,:) = xxxx;
    E_eta2(k,1) = eta_2;
    E_eta2(k,2) = pe;
end

xx = zeros(interval_of_eta_2,1);
xxtemp = E_eta2(:,1); 
for iiiii = 1:length(E_eta2)
    xx(iiiii) = norm(Rot(phi_1)*s1 + Rot(xxtemp(iiiii))*s2 - Rot(phi_1)*t1 - Rot(xxtemp(iiiii))*t2);
end
yy = E_eta2(:,2);
Tosave = [xx,yy];

figure(5);
plot(xx,yy,'color',[1 0 0 1],LineWidth=2.5);
% xlim([0.9 1.25])
% ylim([-0.00001 0.0031])

xlabel('$\lambda(\xi)$','Interpreter','latex','fontsize',18);
ylabel('$E_{act}\left( \xi \right)$','Interpreter','latex','fontsize',18);hold on;
% set ( gca, 'xdir', 'reverse' )

figure(6);
plot((E_eta2(:,1)-phi_1),E_eta2(:,2),'r',LineWidth=2.5);
xlabel('$\xi$','Interpreter','latex','fontsize',18);
ylabel('$E_{act}\left( \xi \right)$','Interpreter','latex','fontsize',18);
% xlim([-0.2 1.2])
% ylim([-0.00001 0.0031])


Tosave_2 = [(E_eta2(:,1)-phi_1),E_eta2(:,2)];

toc