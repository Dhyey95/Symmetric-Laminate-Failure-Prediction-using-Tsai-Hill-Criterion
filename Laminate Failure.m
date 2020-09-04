clear all
clc

th = 0.005 ; %thickness of each ply
% Fiber Properties Input
prompt = 'Material Name: ';
f = input(prompt, 's')
E1_f = 34*10^6;
E2_f = 2.2*10^6;
mu_f_12 = 0.2;
mu_f_21 = mu_f_12*E2_f/E1_f;
G12_f = 4*10^6;
cf = 0.54;

%Matrix Properties Input
prompt = 'Material Name: ';
m = input(prompt, 's')
Em = 0.62*10^6;
mu_m = 0.35;
cm = 0.46;
Gm = 0.24*10^6;

% Hybrid VR Model
E1_H = cf*E1_f + cm*Em;
a1 = ((cf/E2_f) + (cm/Em));
a2 = (cf*cm)*(((mu_f_12*Em)-(mu_m*E1_f))^2)/((E1_f*Em)*(cm*Em + cf*E1_f));
E2_H = 1/(a1-a2);
G12_H = 1/((cf/G12_f)+(cm/Gm));
v12_H = cf*mu_f_12 + cm*mu_m;
v21_H = v12_H*E2_H/E1_H;
q_11 = (E1_H/(1-v12_H*v21_H));
q_12 = ((v21_H*E1_H)/(1-v12_H*v21_H));
q_22 = (E2_H/(1-v12_H*v21_H));
q_66 = G12_H;
Qvr = [q_11 q_12 0; q_12 q_22 0; 0 0 q_66];
E1vr = ((q_11*q_22)-(q_12^2))/q_22;
E2vr = ((q_11*q_22)-(q_12^2))/q_11;
v12vr = q_12/q_22;
v21vr = v12vr*E2vr/E1vr;
G12vr = q_66;
E1 = E1vr;
E2 = E2vr;
G12 = G12vr;
v12 = v12vr;
v21 = v21vr;
elastic_prop = {'E1', 'E2', 'G12', 'v12', 'v21'};
properties_table = table(E1,E2,G12,v12,v21,'VariableNames',elastic_prop)
disp(properties_table)

ce_f_1 = -0.3*10^-6; % in inverse degree farenheit
ce_f_2 = 8.3*10^-6 ; % in inverse degree farenheit
ce_m = 25*10^-6; % in inverse degree farenheit
ce_1 = ((E1_f*ce_f_1*cf)+(Em*ce_m*cm))/(E1_f*cf + Em*cm);
ce_2 = ce_f_1*cf*(1+mu_f_12) + ce_m*cm*(1+mu_m) - v12*ce_1;
a = [ce_1; ce_2; 0];
c_T = 300;
T_T = 72;
dT = T_T - c_T;

% Strength Paramaters
Ft = 535*10^3;
Fmt = 10*10^3;
Fmc = 30*10^3;
Fms = 15*10^3;
Ks = (1-cf*(1-Em/E2_f))/(1-(4*cf/pi())^0.5*(1-(Em/E2_f)));
Kt = (1-cf*(1-Gm/G12_f))/(1-(4*cf/pi())^0.5*(1-(Gm/G12_f)));
F1_t = Ft*(cf + cm*(Em/E2_f));
F1_c = Gm/(1-cf);
F2_c = Fmc/Ks;
F2_t = Fmt/Ks;
F6  = Fms/Kt;
X = F1_t; 
Xc = -F1_c; 
Y = F2_t; 
Yc = -F2_c; 
S = F6;


% stack up 
n = double(input('Layers: '));
theta = zeros(n,1);
for i = 1:n
    theta(i) = double(input('fiber orientation: '));
end
z = zeros(n,1);
z1 = zeros(3,3);
z2 = zeros(3,1);
Qa = cat(3,z1);
Ta = cat(3,z1);
Teg = cat(3,z1);
Nt = zeros(3,1);
Mt = zeros(3,1);
al = cat(3,z2);
t_strain = cat(3,z2);
tpg = zeros(n,1);
btg = zeros(n,1);
d_load = 10;
u_load = 10000
load = 0:d_load:u_load;
m_strain = zeros(3,size(load,2))
m_K = zeros(3,size(load,2))
mark_fail = zeros(n,1);
A_1 = zeros(3,3);
B_1 = zeros(3,3);
D_1 = zeros(3,3);
number = 1;
failed_ply = zeros(n,1);
load_failure = zeros(n,2);
lim_eq = cat(n,1)
for j =  0:d_load:u_load
    N1 = load(1,number);
    Nm = [N1;0;0];
    Mm = [0;0;0];    
    A_1 = zeros(3,3);
    B_1 = zeros(3,3);
    D_1 = zeros(3,3);
    for i = 1:n
        tp = -0.03 + (i-1)*th ;
        tpg(i) = tp;
        bt = -0.035 + (i-1)*th;
        btg(i) = bt;
        z(i) = (tp+bt)/2; 
        t = theta(i);
        c = cosd(t);
        s = sind(t);
        if mark_fail(i,1) == i
            E1 = .00001;
            E2 = .00001;
            G12 = .00001;
        else
            E1 = E1vr;
            E2 = E2vr;
            G12 = G12vr;
        end
        q11 = E1/(1 - (v12)*(v21));
        q12 = v12*E2/(1 - (v12)*(v21));
        q21 = v12*E2/(1 - (v12)*(v21));
        q22 = E2/(1-(v12)*(v21));
        q66 = G12;
        q_local = [q11 q12 0; q21 q22 0; 0 0 q66];
        Ts= [c*c s*s 2*c*s ; s*s c*c -2*c*s ; -c*s  c*s  (c*c)-(s*s)];
        Ta(:,:,i) = Ts;
        Te = [c^2 s^2 c*s ; s^2 c^2 -c*s ; -2*c*s  2*c*s  (c*c)-(s*s)];
        Teg(:,:,i) = Te;
        Q  = inv(Ts)*q_local*Te; 
        Qa(:,:,i) = Q ;
        A_1 = A_1 + Q*th;
        B_1 = B_1 + (0.5).*Q*(tp*tp - bt*bt); %[B] matrix
        D_1 = D_1 + (1/3).*Q*(tp*tp*tp - bt*bt*bt); %[D] matrix
        al = Ts'*a*dT;
        al_a(:,1,i) = al;
        Nt = Nt + Q*al*(tp - bt);
        Mt = Mt + (1/2)*Q*al*(tp*tp - bt*bt);
        t_strain(:,1,i) = al_a(:,1,i);
        Sc = [1/E1 -(v12)/E1 0; -(v12)/E1 1/E2 0; 0 0 1/G12];
    end
    A = A_1;
    B = B_1;
    D = D_1;
    
    if j == d_load
        variable = {'A_Matrix', 'B_Matrix', 'D_Matrix'}
        abd_mat = table(A,B,D,'VariableNames',variable)
        disp(abd_mat)
    end
    
    C = zeros(6,6);
    C(1:3,1:3) = A;
    C(1:3,4:6) = B;
    C(4:6,1:3) = B;
    C(4:6,4:6)= D;
    C_inv = inv(C);
    alpha = C_inv(1:3,1:3);
    beta = C_inv(1:3, 4:6);
    delta = C_inv(4:6,4:6);
    
    e_t = alpha*Nt + beta*Mt;
    if j == d_load
        variable = {'thermal_load','residual_strain'};
        t_str_strain = table(Nt,e_t,'VariableNames', variable);
        disp(t_str_strain)
    end
    k_t = (beta)'*Nt + delta*Mt;
    m_strain(:,number) = alpha*Nm + beta*Mm ;
    m_K(:,number) = (beta)'*Nm + delta*Mm;
    g_sigma = cat(3,z2);
    l_sigma = cat(3,z2);
    t_sigma = cat(3,z2);
    l_strain = cat(3,z2);
    total_sigma = cat(3,z2);
    for r = 1:n
        g_sigma(:,1,r) = Qa(:,:,r)*(m_strain(:,number) + z(r).*m_K(:,number));
        t_sigma(:,1,r) = Qa(:,:,r)*(e_t-t_strain(:,1,r)); %Moment effect is ignored
        
        if j == d_load
            variable = {strcat('Thermal_stress_in_ply',sprintf('%d',r)),};
            thermal_stresses = table(t_sigma(:,1,r),'VariableNames',variable);
            disp(thermal_stresses)
        end

        total_sigma(:,1,r) = g_sigma(:,1,r);
        l_sigma(:,1,r) = Ta(:,:,r)*total_sigma(:,1,r);
        l_strain(:,1, r) = Sc*l_sigma(:,1,r);
        lim_eq(r,1) = l_sigma(1,1,r)^2/X^2 + l_sigma(2,1,r)^2/Y^2 + l_sigma(3,1,r)^2/S^2 - (l_sigma(1,1,r)*l_sigma(2,1,r))/X^2;
        if lim_eq(r,1) > 1
            mark = r;
            fail = theta(mark);
            Nu = Nm(1);
            failed_ply(mark,1) = theta(mark);
            mark_fail(mark,1) = mark;
            load_failure(mark,1) = mark;
            load_failure(mark,2) = Nm(1);
            pp = 'failure of ply number %d with fiber orientation %d degree at load %.f psi \n';
            fprintf(pp, mark, fail, Nu)
        end
    end
    
    number = number+1;
    if sum(mark_fail) == size(mark_fail,1)*(mark_fail(1,1) + mark_fail(n,1))/2 && sum(mark_fail) ~= 0 
        disp('all the plies have failed')
        break
    end
end

plot(m_strain(1,:),load, 'r', 'DisplayName', 'Stress-Strain Curve')
xlabel('mid-plane strain in x-direction')
ylabel('Nx in lb/inch')
title('stress vs strain curve')

