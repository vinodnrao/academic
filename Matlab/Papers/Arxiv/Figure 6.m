% Matlab Code for Figure 6 in `Characterisation of a satellite-to-ground channel for continuous variable quantum key distribution protocol'

clear all;

function gx=h(x)
    if x==1
        gx=0;
    else
    gx=((x+1)/2).*log2((x+1)/2)-((x-1)/2).*log2((x-1)/2);
    end;
end;

% Bypass Channel Key Rate

mu=300;
eta_s=0; %eta here represents eta_T in the paper
c=sqrt(mu.*mu-1);
eta_AE=.05;
eta_t = .1;%1.5*eta_s.*(1-eta_AE)
ex = 0.001;
beta = 1;
eta_d = 1;

k=0;
for Teq = 0:0.00000001:0.0039
    k=k+1;
    
    T=((sqrt(Teq)-sqrt((1-eta_t).*eta_s.*(1-eta_AE))).^2./eta_t./eta_AE.*eta_d); %T represents eta_E in the text
    mu_E=ex./eta_t./(1-T) +1;
    c_E = sqrt(mu_E.*mu_E-1);
    
    CM =[mu,0,sqrt(eta_t).* sqrt(T).* sqrt(eta_AE).* c + sqrt(1-eta_t) .* sqrt(eta_s) .* sqrt(1-eta_AE) .*c,0,0,0,-sqrt(1-T) .* sqrt(eta_AE).* c,0;
        0,mu,0,-sqrt(eta_t) .* sqrt(T).* sqrt(eta_AE).* c- sqrt(1-eta_t).* sqrt(eta_s).* sqrt(1-eta_AE).* c,0,0,0,sqrt(1-T) .* sqrt(eta_AE) .* c;
        sqrt(eta_t) .* sqrt(T) .* sqrt(eta_AE).* c + sqrt(1-eta_t).* sqrt(eta_s).* sqrt(1-eta_AE).* c, 0,(sqrt(eta_t).* (T.* (eta_AE .*mu-eta_AE+1)+(1-T) .*mu_E)-sqrt(1-eta_t) .*sqrt(eta_s) .*sqrt(T) .*(-sqrt(eta_AE) .*mu.* sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE))) .*sqrt(eta_t)-(sqrt(eta_t) .*sqrt(eta_s) .*sqrt(T) .*(-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE).* sqrt(eta_AE))-sqrt(1-eta_t) .*(1+eta_s .*((1-eta_AE) .*mu+eta_AE)-eta_s)) .*sqrt(1-eta_t),0,sqrt(eta_t) .*sqrt(1-T) .*c_E,0,sqrt(eta_t) .*(-sqrt(T) .*(eta_AE .*mu-eta_AE+1) .*sqrt(1-T)+sqrt(1-T).* mu_E .*sqrt(T))+sqrt(1-eta_t) .*sqrt(eta_s) .*sqrt(1-T) .*(-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE)),0;
        0,-sqrt(eta_t) .*sqrt(T) .*sqrt(eta_AE).* c-sqrt(1-eta_t).* sqrt(eta_s).* sqrt(1-eta_AE).* c,0,(sqrt(eta_t) .*(T .*(eta_AE .*mu-eta_AE+1)+(1-T) .*mu_E)-sqrt(1-eta_t) .*sqrt(eta_s).* sqrt(T) .*(-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE))) .*sqrt(eta_t)-(sqrt(eta_t).* sqrt(eta_s) .*sqrt(T).* (-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE))-sqrt(1-eta_t) .*(1+eta_s .*((1-eta_AE) .*mu+eta_AE)-eta_s)) .*sqrt(1-eta_t),0,-sqrt(eta_t) .*sqrt(1-T) .*c_E,0,sqrt(eta_t) .*(-sqrt(T) .*(eta_AE .*mu-eta_AE+1) .*sqrt(1-T)+sqrt(1-T) .*mu_E .*sqrt(T))+sqrt(1-eta_t) .*sqrt(eta_s) .*sqrt(1-T) .*(-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE));
        0,0,sqrt(eta_t) .*sqrt(1-T).*c_E,0,mu_E,0,sqrt(T) .*c_E,0;
        0,0,0,-sqrt(eta_t) .*sqrt(1-T).*c_E,0,mu_E,0,-sqrt(T).*c_E;
        -sqrt(1-T) .*sqrt(eta_AE).* c,0,sqrt(eta_t) .*(-sqrt(T) .*(eta_AE.* mu-eta_AE+1).* sqrt(1-T)+sqrt(1-T) .*mu_E .*sqrt(T))+sqrt(1-eta_t) .*sqrt(eta_s) .*sqrt(1-T) .*(-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE)),0,sqrt(T) .*c_E,0,(1-T) .*(eta_AE .*mu-eta_AE+1)+T .*mu_E,0;
        0,sqrt(1-T).* sqrt(eta_AE).* c,0,sqrt(eta_t).* (-sqrt(T) .*(eta_AE.* mu-eta_AE+1).* sqrt(1-T)+sqrt(1-T) .*mu_E .*sqrt(T))+sqrt(1-eta_t) .*sqrt(eta_s) .*sqrt(1-T).* (-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE)),0,-sqrt(T) .*c_E,0,(1-T) .*(eta_AE .*mu-eta_AE+1)+T .*mu_E];
    
    CM_E = CM(5:8,5:8);
    sigma_EB = CM(5:8,3:4);
    CM_ECB = CM_E - 1/CM(3,3)* sigma_EB* [1 0;0 0]* transpose(sigma_EB);
    
    OM = [0 1.0 0 0 0 0 0 0;
        -1 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0;
        0 0 -1 0 0 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 -1 0 0 0;
        0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 -1 0];
    
    e_E = sort(imag(eig(OM(1:4,1:4)*CM_E)));
    e_ECB = sort(imag(eig(OM(1:4,1:4)*CM_ECB)));
    
    v_E = e_E(length(e_E)/2+1:length(e_E));
    v_ECB = e_ECB(length(e_ECB)/2+1:length(e_ECB));
    
    for i=1:2
        if v_E(i)<1;
            1-v_E(i);
            v_E(i)=1;
        end;
        if v_ECB(i)<1;
            1-v_ECB(i);
            v_ECB(i)=1;
        end;
    end;
    
    X_EB(k) = h(v_E(1))+h(v_E(2))-h(v_ECB(1))-h(v_ECB(2));
    
    CM_A = CM(1:2,1:2);
    CM_AB = CM(1:4,1:4);
    sigma_AB = CM(1:2,3:4);
    CM_ACB = CM_A - 1/CM(3,3)* sigma_AB* [1 0;0 0]* transpose(sigma_AB);
    I_AB(k) = (1/2)*log2((CM_A(1,1)+1)*CM(3,3)/((CM_A(1,1)+1)*CM(3,3)-CM(1,3)*CM(3,1)));
   
    R(k) = beta*I_AB(k) - X_EB(k);

end;

figure(1)
semilogy(0:0.00000001:0.0039,R,'LineWidth',1.5,'Color','b')
xticks([0.00031 0.00039 0.0005 0.00063 0.00079 0.001])
xticklabels({'35','34','33','32','31','30'})
set(gca,'Xdir','reverse')
set(gca,'FontSize',14)
xlim([0.00028, 0.001])
ylim([0.000001, 0.01])
xlabel('T_{eq} (/dB)','FontSize',14)
ylabel('Rate (bits per use)','FontSize',14)

hold on

mu=150;
eta_s=0; %eta here represents eta_T in the paper
c=sqrt(mu.*mu-1);
eta_AE=.01;
eta_t = .1;
ex = 0.005;
beta = 0.95;
eta_d = 0.8;

k=0;
for Teq = 0:0.00000001:0.0039
    k=k+1;
    
    T=((sqrt(Teq)-sqrt((1-eta_t).*eta_s.*(1-eta_AE))).^2./eta_t./eta_AE.*eta_d); %T represents eta_E in the text
    mu_E=ex./eta_t./(1-T) +1;
    c_E = sqrt(mu_E.*mu_E-1);
    
    CM =[mu,0,sqrt(eta_t).* sqrt(T).* sqrt(eta_AE).* c + sqrt(1-eta_t) .* sqrt(eta_s) .* sqrt(1-eta_AE) .*c,0,0,0,-sqrt(1-T) .* sqrt(eta_AE).* c,0;
        0,mu,0,-sqrt(eta_t) .* sqrt(T).* sqrt(eta_AE).* c- sqrt(1-eta_t).* sqrt(eta_s).* sqrt(1-eta_AE).* c,0,0,0,sqrt(1-T) .* sqrt(eta_AE) .* c;
        sqrt(eta_t) .* sqrt(T) .* sqrt(eta_AE).* c + sqrt(1-eta_t).* sqrt(eta_s).* sqrt(1-eta_AE).* c, 0,(sqrt(eta_t).* (T.* (eta_AE .*mu-eta_AE+1)+(1-T) .*mu_E)-sqrt(1-eta_t) .*sqrt(eta_s) .*sqrt(T) .*(-sqrt(eta_AE) .*mu.* sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE))) .*sqrt(eta_t)-(sqrt(eta_t) .*sqrt(eta_s) .*sqrt(T) .*(-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE).* sqrt(eta_AE))-sqrt(1-eta_t) .*(1+eta_s .*((1-eta_AE) .*mu+eta_AE)-eta_s)) .*sqrt(1-eta_t),0,sqrt(eta_t) .*sqrt(1-T) .*c_E,0,sqrt(eta_t) .*(-sqrt(T) .*(eta_AE .*mu-eta_AE+1) .*sqrt(1-T)+sqrt(1-T).* mu_E .*sqrt(T))+sqrt(1-eta_t) .*sqrt(eta_s) .*sqrt(1-T) .*(-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE)),0;
        0,-sqrt(eta_t) .*sqrt(T) .*sqrt(eta_AE).* c-sqrt(1-eta_t).* sqrt(eta_s).* sqrt(1-eta_AE).* c,0,(sqrt(eta_t) .*(T .*(eta_AE .*mu-eta_AE+1)+(1-T) .*mu_E)-sqrt(1-eta_t) .*sqrt(eta_s).* sqrt(T) .*(-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE))) .*sqrt(eta_t)-(sqrt(eta_t).* sqrt(eta_s) .*sqrt(T).* (-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE))-sqrt(1-eta_t) .*(1+eta_s .*((1-eta_AE) .*mu+eta_AE)-eta_s)) .*sqrt(1-eta_t),0,-sqrt(eta_t) .*sqrt(1-T) .*c_E,0,sqrt(eta_t) .*(-sqrt(T) .*(eta_AE .*mu-eta_AE+1) .*sqrt(1-T)+sqrt(1-T) .*mu_E .*sqrt(T))+sqrt(1-eta_t) .*sqrt(eta_s) .*sqrt(1-T) .*(-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE));
        0,0,sqrt(eta_t) .*sqrt(1-T).*c_E,0,mu_E,0,sqrt(T) .*c_E,0;
        0,0,0,-sqrt(eta_t) .*sqrt(1-T).*c_E,0,mu_E,0,-sqrt(T).*c_E;
        -sqrt(1-T) .*sqrt(eta_AE).* c,0,sqrt(eta_t) .*(-sqrt(T) .*(eta_AE.* mu-eta_AE+1).* sqrt(1-T)+sqrt(1-T) .*mu_E .*sqrt(T))+sqrt(1-eta_t) .*sqrt(eta_s) .*sqrt(1-T) .*(-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE)),0,sqrt(T) .*c_E,0,(1-T) .*(eta_AE .*mu-eta_AE+1)+T .*mu_E,0;
        0,sqrt(1-T).* sqrt(eta_AE).* c,0,sqrt(eta_t).* (-sqrt(T) .*(eta_AE.* mu-eta_AE+1).* sqrt(1-T)+sqrt(1-T) .*mu_E .*sqrt(T))+sqrt(1-eta_t) .*sqrt(eta_s) .*sqrt(1-T).* (-sqrt(eta_AE) .*mu .*sqrt(1-eta_AE)+sqrt(1-eta_AE) .*sqrt(eta_AE)),0,-sqrt(T) .*c_E,0,(1-T) .*(eta_AE .*mu-eta_AE+1)+T .*mu_E];
    
    CM_E = CM(5:8,5:8);
    sigma_EB = CM(5:8,3:4);
    CM_ECB = CM_E - 1/CM(3,3)* sigma_EB* [1 0;0 0]* transpose(sigma_EB);
    
    OM = [0 1.0 0 0 0 0 0 0;
        -1 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0;
        0 0 -1 0 0 0 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 -1 0 0 0;
        0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 -1 0];
    
    e_E = sort(imag(eig(OM(1:4,1:4)*CM_E)));
    e_ECB = sort(imag(eig(OM(1:4,1:4)*CM_ECB)));
    
    v_E = e_E(length(e_E)/2+1:length(e_E));
    v_ECB = e_ECB(length(e_ECB)/2+1:length(e_ECB));
    
    for i=1:2
        if v_E(i)<1;
            1-v_E(i);
            v_E(i)=1;
        end;
        if v_ECB(i)<1;
            1-v_ECB(i);
            v_ECB(i)=1;
        end;
    end;
    
    X_EB(k) = h(v_E(1))+h(v_E(2))-h(v_ECB(1))-h(v_ECB(2));
    
    CM_A = CM(1:2,1:2);
    CM_AB = CM(1:4,1:4);
    sigma_AB = CM(1:2,3:4);
    CM_ACB = CM_A - 1/CM(3,3)* sigma_AB* [1 0;0 0]* transpose(sigma_AB);
    I_AB(k) = (1/2)*log2((CM_A(1,1)+1)*CM(3,3)/((CM_A(1,1)+1)*CM(3,3)-CM(1,3)*CM(3,1)));
   
    R(k) = beta*I_AB(k) - X_EB(k);

end;

figure(1)
semilogy(0:0.00000001:0.0039,R,'LineWidth',1.5,'Color','r')
set(gca,'Xdir','reverse')
set(gca,'FontSize',14)

qw{1} = plot(nan, 'b');
qw{2} = plot(nan, 'r');
qw{3} = plot(nan, 'k--');
legend([qw{:}], {'Ideal Res.','Realistic Res.','Unrestricted'}, 'location', 'northwest')

hold on

v = 2.4;
a = v;
be = 0.95;
e = 1e-4;

len = 0:0.01:200;
% Define t(l) as a function handle
%t = 1 ./ 10.^(0.02 * l);
tr = @(len) 1 ./ 10.^(0.02 * len);

% Define b(l)
b = @(len) (tr(len) .* (v - 1) + 1 + e);

% Define c(l)
c = @(len) sqrt(tr(len) .* (v^2 - 1));

% Define z(l)
z = @(len) sqrt((a + b(len)).^2 - 4 .* c(len).^2);

% Define f(E)
f = @(E) ((E + 1)/2) .* log2((E + 1)/2) - ((E - 1)/2) .* log2((E - 1)/2);

% Define E1, E2, E3
E1 = @(len) 0.5 * (z(len) + (b(len) - a));
E2 = @(len) 0.5 * (z(len) - (b(len) - a));
E3 = @(len) sqrt(a * (a - (c(len).^2 ./ b(len))));

% Define fE1, fE2, fE3, fE12
fE1 = @(len) f(E1(len));
fE2 = @(len) f(E2(len));
fE3 = @(len) f(E3(len));
fE12 = @(len) fE1(len) + fE2(len);

% Define Chi1
Chi1 = @(len) fE12(len) - fE3(len);

% Define SNR and IAB
SNR = @(len) (tr(len) .* (v - 1)) ./ (1 + e);
IAB = @(len) 0.5 * log2(1 + SNR(len));

% Define r(l) and r2(l)
r = @(len) (be * IAB(len) - Chi1(len));
r2 = @(len) 2e6 * r(len);

axes('Position',[.2 .2 .25 .25])
box on
semilogy(tr(len), r(len), "black", 'LineWidth', 1, 'LineStyle','--');
xticks([0.001 0.00316 0.01])
xticklabels({'30','25','20'})
set(gca,'Xdir','reverse');
xlim([0.001, 0.01]);
%ylim([0.001, 0.05]);

hold off

hold off

saveas(figure(1),'Keyrate_Emma_2b.pdf')
print(gcf,'Keyrate_Emma_2b.jpeg','-djpeg','-r1000'); 
