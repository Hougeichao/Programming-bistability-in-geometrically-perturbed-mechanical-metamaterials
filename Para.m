function[s1,s2,s3,s4,t1,t2,t3,t4,u1,u2,u3,u4,v1,v2,v3,v4,Module_1_xi,Module_1_lambda,Module_2_xi,Module_2_lambda,Coor_ref,Coor_def] = Para(s11,s12,s21,s22,phi_1,phi_2,phi_3,phi_4,l1R,l2R,l1D,l2D,p)

s1 = [s11;s12];s2 = [s21;s22];Ident = Rot(0);
R1 = Rot(phi_1);R2 = Rot(phi_2); R3 = Rot(phi_3); R4 = Rot(phi_4);
delta_12 = R1 - R2;
delta_23 = R2 - R3;
delta_34 = R3 - R4;
delta_41 = R4 - R1;

%% matrix that solves all the equality constraints
s3 = (delta_34)\delta_41*s1 + (-Ident - (delta_34)\delta_23)*s2;
s4 = (-Ident - (delta_34)\delta_41)*s1 + (delta_34)\delta_23*s2;
t1 = Ident*s1 + (delta_12)\R2*l1R - (delta_12)\l1D;
t2 = Ident*s2 - (delta_12)\R1*l1R + (delta_12)\l1D;
t3 = (delta_34)\delta_41*s1 - (Ident + (delta_34)\delta_23)*s2 - (delta_34)\R4*l1R + (delta_34)\l1D;
t4 = (-Ident - (delta_34)\delta_41)*s1 + (delta_34)\delta_23*s2 + (delta_34)\R3*l1R - (delta_34)\l1D;
u1 = Ident*s1 + (delta_41)\R4*l2R - (delta_41)\l2D;
u2 = Ident*s2 + (delta_23)\R3*l2R - (delta_23)\l2D;
u3 = (delta_34)\delta_41*s1 + (-Ident - (delta_34)\delta_23)*s2 - (delta_23)\R2*l2R + (delta_23)\l2D;
u4 = -(Ident + (delta_34)\delta_41)*s1 + (delta_34)\delta_23*s2 - (delta_41)\R1*l2R + (delta_41)\l2D;
v1 = Ident*s1 + (delta_41)\R4*l2R + (delta_12)\R2*l1R - (delta_41)\l2D - (delta_12)\l1D;
v2 = Ident*s2 + (delta_23)\R3*l2R - (delta_12)\R1*l1R - (delta_23)\l2D + (delta_12)\l1D;
v3 = (delta_34)\delta_41*s1 - (Ident + (delta_34)\delta_23)*s2 - (delta_23)\R2*l2R - (delta_34)\R4*l1R + (delta_23)\l2D + (delta_34)\l1D;
v4 = -(Ident + (delta_34)\delta_41)*s1 + (delta_34)\delta_23*s2 - (delta_41)\R1*l2R + (delta_34)\R3*l1R + (delta_41)\l2D - (delta_34)\l1D;

%% Module calculation
S = [-1/2*eye(2),3/8*eye(2),1/8*eye(2),-1/4*eye(2);1/4*eye(2),-1/8*eye(2),-3/8*eye(2),1/2*eye(2);-1/8*eye(2),1/2*eye(2),1/4*eye(2),-3/8*eye(2);3/8*eye(2),-1/4*eye(2),-1/2*eye(2),1/8*eye(2)];
G = 2*(S')*S;
k1_1 = 2*dot([s2;t2;u2;v2],G*[s2;t2;u2;v2]);
k1_vector_1 = 2*[s3' t3' u3' v3';s4' t4' u4' v4']*G*[s2;t2;u2;v2];
K1_tensor_1 = 2*[s3' t3' u3' v3';s4' t4' u4' v4']*G*[s3 s4;t3 t4;u3 u4;v3 v4];
k1_2 = 2*dot([R2*s2;R2*t2;R2*u2;R2*v2],G*[R2*s2;R2*t2;R2*u2;R2*v2]);
k1_vector_2 = 2*[(R3*s3)' (R3*t3)' (R3*u3)' (R3*v3)';(R4*s4)' (R4*t4)' (R4*u4)' (R4*v4)']*G*[R2*s2;R2*t2;R2*u2;R2*v2];
K1_tensor_2 = 2*[(R3*s3)' (R3*t3)' (R3*u3)' (R3*v3)';(R4*s4)' (R4*t4)' (R4*u4)' (R4*v4)']*G*[R3*s3 R4*s4;R3*t3 R4*t4;R3*u3 R4*u4;R3*v3 R4*v4];
C1 = (norm(s1 - t1))^2 + norm(s2 - t2)^2; C2 = dot((s1 - t1),(s2 - t2));C3 = dot((s1 - t1),Rot(pi/2)*(s2 - t2));

Module_1_xi = k1_1 - dot(k1_vector_1,inv(K1_tensor_1)*k1_vector_1);
Module_1_lambda = Module_1_xi/(C3/sqrt(C1+2*C2))^2;
Module_2_xi = k1_2 - dot(k1_vector_2,inv(K1_tensor_2)*k1_vector_2);
Module_2_lambda = Module_2_xi/((cos(phi_2 - phi_1)*C3 - sin(phi_2 - phi_1)*C2)/sqrt(C1 + 2*sin(phi_2 - phi_1)*C3 + 2*cos(phi_2 - phi_1)*C2))^2;

%% Coordinate storage
xr_11 = [0;0]; xr_21 = s1; xr_31 = s1 - u1; xr_41 = t1;
xr_12 = s1; xr_22 = s1 + s2; xr_32 = s1 + s2 - t2; xr_42 = s1 + u2;
xr_13 = s1 + s2; xr_23 = s1 + s2 + t3; xr_33 = s1 + s2 + t3 - v3; xr_43 = -s4;
xr_14 = [0;0]; xr_24 = -s4; xr_34 = -s4 + u4; xr_44 = -t4;

yr_11 = [0;0]; yr_21 = R1*s1; yr_31 = R1*(s1 - u1); yr_41 = R1*t1;
yr_12 = R1*s1; yr_22 = R1*s1 + R2*s2; yr_32 = R1*s1 + R2*(s2 - t2); yr_42 = R1*s1 + R2*u2;
yr_13 = R1*s1 + R2*s2; yr_23 = R1*s1 + R2*s2 + R3*t3; yr_33 = R1*s1 + R2*s2 + R3*(t3 - v3); y_43 = -R4*s4;
yr_14 = [0;0]; yr_24 = -R4*s4; yr_34 = R4*(-s4 + u4); yr_44 = -R4*t4;

x_r = [xr_11(1) xr_12(1) xr_13(1) xr_14(1);xr_21(1) xr_22(1) xr_23(1) xr_24(1);xr_31(1) xr_32(1) xr_33(1) xr_34(1);xr_41(1) xr_42(1) xr_43(1) xr_44(1)];
y_r = [xr_11(2) xr_12(2) xr_13(2) xr_14(2);xr_21(2) xr_22(2) xr_23(2) xr_24(2);xr_31(2) xr_32(2) xr_33(2) xr_34(2);xr_41(2) xr_42(2) xr_43(2) xr_44(2)];
Coor_ref = [xr_11';xr_41';xr_31';xr_21';xr_42';xr_32';xr_22';xr_23';xr_33';xr_43';xr_34';xr_44'];

x_d = [yr_11(1) yr_12(1) yr_13(1) yr_14(1);yr_21(1) yr_22(1) yr_23(1) yr_24(1);yr_31(1) yr_32(1) yr_33(1) yr_34(1);yr_41(1) yr_42(1) y_43(1) yr_44(1)];
y_d = [yr_11(2) yr_12(2) yr_13(2) yr_14(2);yr_21(2) yr_22(2) yr_23(2) yr_24(2);yr_31(2) yr_32(2) yr_33(2) yr_34(2);yr_41(2) yr_42(2) y_43(2) yr_44(2)];
Coor_def = [yr_11';yr_41';yr_31';yr_21';yr_42';yr_32';yr_22';yr_23';yr_33';y_43';yr_34';yr_44'];

%% Plot unit cell and 2-by-2 patterns
if p == 1
    xmin_unit = min(min([x_r;x_d]));xmax_unit = max(max([x_r;x_d]));
    ymin_unit = min(min([y_r;y_d]));ymax_unit = max(max([y_r;y_d]));
    figure(1);hold on;
    xlim([xmin_unit xmax_unit])
    ylim([ymin_unit ymax_unit])
    patch(x_r,y_r,[200/255 200/255 200/255],'FaceAlpha',0.78,'LineWidth',1.2);axis equal;axis off;
    figure(2);hold on;
    xlim([xmin_unit xmax_unit])
    ylim([ymin_unit ymax_unit])
    patch(x_d,y_d,[200/255 200/255 200/255],'FaceAlpha',0.78,'LineWidth',1.2);axis equal;axis off;
    
    % plot 2-by-2 pattern
    x_r_2 = [xr_11(1) xr_12(1) xr_13(1) xr_14(1);xr_21(1) xr_22(1) xr_23(1) xr_24(1);xr_31(1) xr_32(1) xr_33(1) xr_34(1);xr_41(1) xr_42(1) xr_43(1) xr_44(1)]+l1R(1);
    y_r_2 = [xr_11(2) xr_12(2) xr_13(2) xr_14(2);xr_21(2) xr_22(2) xr_23(2) xr_24(2);xr_31(2) xr_32(2) xr_33(2) xr_34(2);xr_41(2) xr_42(2) xr_43(2) xr_44(2)]+l1R(2);
    x_r_3 = [xr_11(1) xr_12(1) xr_13(1) xr_14(1);xr_21(1) xr_22(1) xr_23(1) xr_24(1);xr_31(1) xr_32(1) xr_33(1) xr_34(1);xr_41(1) xr_42(1) xr_43(1) xr_44(1)]+l2R(1);
    y_r_3 = [xr_11(2) xr_12(2) xr_13(2) xr_14(2);xr_21(2) xr_22(2) xr_23(2) xr_24(2);xr_31(2) xr_32(2) xr_33(2) xr_34(2);xr_41(2) xr_42(2) xr_43(2) xr_44(2)]+l2R(2);
    x_r_4 = [xr_11(1) xr_12(1) xr_13(1) xr_14(1);xr_21(1) xr_22(1) xr_23(1) xr_24(1);xr_31(1) xr_32(1) xr_33(1) xr_34(1);xr_41(1) xr_42(1) xr_43(1) xr_44(1)]+l1R(1)+l2R(1);
    y_r_4 = [xr_11(2) xr_12(2) xr_13(2) xr_14(2);xr_21(2) xr_22(2) xr_23(2) xr_24(2);xr_31(2) xr_32(2) xr_33(2) xr_34(2);xr_41(2) xr_42(2) xr_43(2) xr_44(2)]+l1R(2)+l2R(2);

    x_d_2 = [yr_11(1) yr_12(1) yr_13(1) yr_14(1);yr_21(1) yr_22(1) yr_23(1) yr_24(1);yr_31(1) yr_32(1) yr_33(1) yr_34(1);yr_41(1) yr_42(1) y_43(1) yr_44(1)]+l1D(1);
    y_d_2 = [yr_11(2) yr_12(2) yr_13(2) yr_14(2);yr_21(2) yr_22(2) yr_23(2) yr_24(2);yr_31(2) yr_32(2) yr_33(2) yr_34(2);yr_41(2) yr_42(2) y_43(2) yr_44(2)]+l1D(2);
    x_d_3 = [yr_11(1) yr_12(1) yr_13(1) yr_14(1);yr_21(1) yr_22(1) yr_23(1) yr_24(1);yr_31(1) yr_32(1) yr_33(1) yr_34(1);yr_41(1) yr_42(1) y_43(1) yr_44(1)]+l2D(1);
    y_d_3 = [yr_11(2) yr_12(2) yr_13(2) yr_14(2);yr_21(2) yr_22(2) yr_23(2) yr_24(2);yr_31(2) yr_32(2) yr_33(2) yr_34(2);yr_41(2) yr_42(2) y_43(2) yr_44(2)]+l2D(2);
    x_d_4 = [yr_11(1) yr_12(1) yr_13(1) yr_14(1);yr_21(1) yr_22(1) yr_23(1) yr_24(1);yr_31(1) yr_32(1) yr_33(1) yr_34(1);yr_41(1) yr_42(1) y_43(1) yr_44(1)]+l1D(1)+l2D(1);
    y_d_4 = [yr_11(2) yr_12(2) yr_13(2) yr_14(2);yr_21(2) yr_22(2) yr_23(2) yr_24(2);yr_31(2) yr_32(2) yr_33(2) yr_34(2);yr_41(2) yr_42(2) y_43(2) yr_44(2)]+l1D(2)+l2D(2);

    xminopt = min(min([x_r;x_r_2;x_r_3;x_r_4;x_d;x_d_2;x_d_3;x_d_4]));
    xmaxopt = max(max([x_r;x_r_2;x_r_3;x_r_4;x_d;x_d_2;x_d_3;x_d_4]));
    yminopt = min(min([y_r;y_r_2;y_r_3;y_r_4;y_d;y_d_2;y_d_3;y_d_4]));
    ymaxopt = max(max([y_r;y_r_2;y_r_3;y_r_4;y_d;y_d_2;y_d_3;y_d_4]));

    figure(3);hold on;axis equal;axis off;
    xlim([xminopt xmaxopt])
    ylim([yminopt ymaxopt])
    patch([x_r,x_r_2,x_r_3,x_r_4],[y_r,y_r_2,y_r_3,y_r_4],[200/255 200/255 200/255],'FaceAlpha',0.78,'LineWidth',1.2);
    
    figure(4); hold on;axis equal;axis off;
    xlim([xminopt xmaxopt])
    ylim([yminopt ymaxopt])
    patch([x_d,x_d_2,x_d_3,x_d_4],[y_d,y_d_2,y_d_3,y_d_4],[200/255 200/255 200/255],'FaceAlpha',0.78,'LineWidth',1.2);
end
end