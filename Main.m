ccc;tic
e1 = [1;0]; e2 = [0;1];

%% define Bravais Lattice vectors
l1R = 1*e1; l2R = 1*e2;
l1D = 1.2*e1; l2D = 0.8*e2;

%% start optimization
options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',10000,'MaxIterations',1000,'ConstraintTolerance',1e-8,'FunctionTolerance',1e-6,'Display','iter');
ic = [0.45;-0.05;0.45;0.05;-0.2*pi;0.15*pi;-0.2*pi;0.2*pi];
[x,y] = fmincon(@(x)Obj(x,l1R,l2R,l1D,l2D),ic,[],[],[],[],[],[],@(x)Constraint(x,l1R,l2R,l1D,l2D),options); 

%% plot figures and post-process
[s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4,Module_1_xi,Module_1_lambda,Module_2_xi,Module_2_lambda,Coor_ref,Coor_def] = Para(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),l1R,l2R,l1D,l2D,1);
[Eb,Ebpp1,Ebpp2,Ebpplambda1,Ebpplambda2] = EnergyBarrier(x,l1R,l2R,l1D,l2D);

toc



