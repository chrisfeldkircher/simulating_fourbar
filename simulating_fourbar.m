%% Declaration of Vars
syms l1_s l2_s l3_s l4_s w1_s w2_s w3_s w4_s k_s
l1 = 50;
l2 = 10;
l3 = 30;
w3 = pi/2;
k = 20;
w1 = (0:0.1:2*pi)';

T_AB = [cos(w1_s) -sin(w1_s) l1_s; sin(w1_s) cos(w1_s) 0; 0 0 1];
T_BC = [cos(w2_s) -sin(w2_s) l2_s; sin(w2_s) cos(w2_s) 0; 0 0 1];
T_CK = [cos(w2_s) -sin(w2_s) k_s; sin(w2_s) cos(w2_s) 0; 0 0 1];
T_CD = [cos(w3_s) -sin(w3_s) l3_s; sin(w3_s) cos(w3_s) 0; 0 0 1];
T_DA = [cos(w4_s) -sin(w4_s) l4_s; sin(w4_s) cos(w4_s) 0; 0 0 1];

%% Solving the System
E=T_AB*T_BC*T_CD*T_DA;

assume(l1_s>0 & l2_s>0 & l3_s>0 & l4_s>0)
assumeAlso(l1_s, 'real')
assumeAlso(l2_s, 'real')
assumeAlso(l3_s, 'real')
assumeAlso(l4_s, 'real')
assumptions;

sol = solve(E(1,3)==0, E(2,3)==0, l4_s, w2_s ,'Real', true);

A = [0;0;1];
B = T_AB*A;
C1 = T_AB*T_BC*A;
K1 = T_AB*T_BC*T_CK*A;
D1 = T_AB*T_BC*T_CD*A;

%% Ddclaring the velocity vars
vel_A = diff(A, w4_s);
vel_B = diff(B, w1_s);
vel_C = diff(C1, w1_s);
vel_K = diff(K1, w1_s);
vel_D_w1 = diff(D1, w1_s);
vel_D_w2 = diff(D1, w2_s);

%% Feeding in the Values
l1_s = l1;
l2_s = l2;
l3_s = l3;
l4 = zeros(length(w1),1);
w1_s = w1;
w2 = zeros(length(w1),1);
w3_s = w3;
w4 = zeros(length(w1),1);
k_s = k;

%% Getting the velocity
vel_C1 = zeros(numel(w1),2);
vel_C = eval(subs(vel_C));
vel_C1(:,1) = vel_C(1:length(w1));
vel_C1(:,2) = vel_C(length(w1)+1:2*length(w1));

vel_K1 = zeros(numel(w1),2);
vel_K = subs(vel_K, 'w2_s', sol.w2_s(1));
vel_K = eval(subs(vel_K));
vel_K1(:,1) = vel_K(1:length(w1));
vel_K1(:,2) = vel_K(length(w1)+1:2*length(w1));

vel_D1 = zeros(numel(w1),2);
vel_D_w1 = subs(vel_D_w1, 'w2_s', sol.w2_s(1));
vel_D_w1 = eval(subs(vel_D_w1));
vel_D1(:,1) = vel_D_w1(1:length(w1));
vel_D1(:,2) = vel_D_w1(length(w1)+1:2*length(w1));

vel_D2 = zeros(numel(w1),2);
vel_D_w2 = subs(vel_D_w2, 'w2_s', sol.w2_s(1));
vel_D_w2 = eval(subs(vel_D_w2));
vel_D2(:,1) = vel_D_w2(1:length(w1));
vel_D2(:,2) = vel_D_w2(length(w1)+1:2*length(w1));


%% Calc the Pos. of K
K_1 = zeros(numel(w1),2);
K1 = subs(K1, 'w2_s', sol.w2_s(1));
K1 = subs(K1);
K_1(:,1) = K1(1:length(w1));
K_1(:,2) = K1(length(w1)+1:2*length(w1));
K = K_1;

%% Calc the Pos. of C
C_1 = zeros(numel(w1),2);
C1 = subs(C1, 'w2_s', sol.w2_s(1));
C1 = subs(C1);
C_1(:,1) = C1(1:length(w1));
C_1(:,2) = C1(length(w1)+1:2*length(w1));
C = C_1;

%% Calc the Pos. of B
B = [l1;0;1];
%% Calc the Pos. of D
D_1 = zeros(numel(w1),2);
D1 = subs(D1, 'w2_s', sol.w2_s(1));
D1 = subs(D1);
D_1(:,1) = D1(1:length(w1));
D_1(:,2) = D1(length(w1)+1:2*length(w1));
D = D_1;

%% Calc length of l4
for i=1:numel(w1)
    l4(i) = sqrt(D(i,1)^2 + D(i,2)^2);
end

%% Calc of w4
for i=1:numel(w1)
    w4(i) = pi-asin(D(i,2)/l4(i));
end

%% Calc of w2
slp_BC = zeros(length(w1),1);
slp_CD = zeros(length(w1),1);

for i=1:numel(w1)
    slp_BC(i) = (C(i,2))/(C(i,1)-l1);
    slp_CD(i) = (D(i,2)-C(i,2))/(D(i,1)-C(i,1));
    
    w2(i) = pi - abs(atan(slp_BC(i)) - atan(slp_CD(i)));
end

%% Plotting
figure
axis equal
axis([-2*l4(1) l1+l2+10 -max(l2,l4(1)) max(l2,l4(1))+10])

% v = VideoWriter('sim.avi'); 
% open(v);

for i=1:numel(w1)
    
%     frame = getframe(gcf);
%     writeVideo(v,frame);
    
%%%%%%%%%%%%%%%%%%%%% Plot 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,1,1);
    hold off
    plot([A(1) B(1) C(i,1) D(i,1) A(1)], [A(2) B(2) C(i,2) D(i,2) A(2)])
    hold on
    title('Animation of the Configuration') 
    
    plot(A(1), A(2), 'r*')
    plot(B(1), B(2), 'r*')
    plot(C(i,1), C(i,2), 'r*')
    plot(K(i,1), K(i,2), '*','color',[.5 .4 .7])
    plot(D(i,1), D(i,2), 'r*')
    
    text(A(1), A(2), 'A');
    text(B(1), B(2), 'B');
    text(C(i,1), C(i,2), 'C');
    text(K(i,1), K(i,2), 'K');
    text(D(i,1), D(i,2), 'D');
    
    axis([-l4(1) l1+l2+10 -max(l2,l4(1)) max(l2,l4(1))+10])
    
%%%%%%%%%%%%%%%%%%%%% Plot 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,1,2);
    title('Path of C, K, D')  
    axis([-l4(1) l1+l2+10 -max(l2,l4(1)) max(l2,l4(1))+10])
%     plot([A(1) B(1) C(i,1) D(i,1) A(1)], [A(2) B(2) C(i,2) D(i,2) A(2)])
    pnt_c = text(C(i,1), C(i,2), 'C');
    %pnt_k = text(K(i,1), K(i,2), 'K');
    pnt_d = text(D(i,1), D(i,2), 'D');
    plot(A(1),  A(2), 'b*')
    hold on
    text(A(1), A(2), 'A');
    text(B(1), B(2), 'B');
    plot(B(1),  B(2), 'b*')
    plot(C(i,1),  C(i,2), 'r.')
    %plot(K(i,1),  K(i,2), '.','color',[.5 .4 .7])
    plot(D(i,1),  D(i,2), 'b.')

    
%%%%%%%%%%%%%%%%%%%%% Plot 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,1,3);
    title('Angular Velocity dD/dw1 & dD/dw2') 
    axis([-l4(1) l1+l2+10 -2*max(l2,l4(1)) 2*max(l2,l4(1))+10])
    plot(vel_D2(i,1), vel_D2(i,2), 'r.')
    hold on
    plot(vel_D1(i,1), vel_D1(i,2), 'b.')
    
    drawnow
    pause(0.035);
    
    delete(pnt_c);
    delete(pnt_d);
    %delete(pnt_k);
end

% close(v);