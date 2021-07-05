clear all;
close all;
tic;

n = 1:200%[1 5:5:200];      % number of STAs
W0 = 1:500;         % OCWmin
m = 0;          % number of backoff stages
M = 1;          % number of RA-RUs 
L_MPDU = 2000;  % Byte
Data_Rate = 6.6667*10^6;%65*10^6; %;  

L_PHY = 40*10^-6;      % preamble length (40 usec)
L_Trigger = 100*10^-6; % length of trigger frame
L_BACK = 68*10^-6;     % length of block ack
SIFS = 16*10^-6;

T = 8*L_MPDU/Data_Rate + 3*SIFS + 3*L_PHY + L_BACK+ L_Trigger;

R_eff = zeros(length(W0), length(n)); 
R_Pc = zeros(length(W0), length(n)); 
R_tau = zeros(length(W0), length(n)); 
R_Th = zeros(length(W0), length(n)); 
for j = 1:length(W0)
    for i = 1:length(n)

            t(i)=fzero(@tau,[0,1],[],n(i),W0(j),m,M);   % tau
            p(i)=1-(1-t(i)/M)^(n(i)-1);              % collision probability
            Ps(i) = t(i)*(1-p(i));                   % Ps
            Es(i) = n(i)*Ps(i);                      % expected number of STAs successful contest at a stage
            eff(i) = Es(i)/M;                        % system efficiency

            Ptr(i) = 1-(1-t(i)/M)^n(i);              % transmission probability  
            Ps_(i) = ( n(i)*(t(i)/M)*(1-(t(i)/M))^(n(i)-1) ) / (1-(1-(t(i)/M))^(n(i))); % success probability
        %     Pi(i) = (1-Ptr(i))^M;                    %    
        %     Th(i) = (M*Ptr(i)*Ps_(i)*8*L_MPDU) / (Pi(i)*T +Ptr(i)*Ps_(i)*T + Ptr(i)*(1-Ps_(i))*T);
            Pw(i) = 1-(1-Ptr(i))^(M-1);
            I(i) = (1-Ptr(i))*(1-Pw(i))*T+(1-Ptr(i))*Pw(i)*T;
            Th(i) = M*(Ptr(i)*Ps_(i)*8*L_MPDU) / (I(i) + Ptr(i)*Ps_(i)*T + Ptr(i)*(1-Ps_(i))*T);
            D(i) = 1 / ( t(i)*(1-(t(i)/M))^(n(i)-1) );
            Ds(i)= 1/(1 - (1-t(i)*(1-p(i)))^n(i) ); 
    end
    R_eff(j,:) = eff';
    R_Pc(j,:) = p';
    R_tau(j,:) = t';
    R_Th(j,:) =  Th';
    eff = [];
end

for i = 1:length(n)
   
    Opt_OCW(i) = find(R_eff(:,i) == max(R_eff(:,i)), 1,'first');
    Opt_OCW2(i) = find(R_eff(:,i) == max(R_eff(:,i)), 1,'Last');
    Opt_EFF1(i) = R_eff(Opt_OCW(i), i);
    Opt_EFF2(i) = R_eff(Opt_OCW2(i), i);
    Opt_Pc1(i) = R_Pc(Opt_OCW(i), i);
    Opt_Pc2(i) = R_Pc(Opt_OCW2(i), i);
    Opt_tau1(i) = R_tau(Opt_OCW(i), i);
    Opt_tau2(i) = R_tau(Opt_OCW2(i), i);
    Opt_Th(i) = R_Th(Opt_OCW2(i), i);
end
Opt_Pc2 = Opt_Pc2';
Opt_tau2 = Opt_tau2';
Opt_Th = Opt_Th'*10^-6; % Mb/s

figure;
hold on;
plot(n, Opt_OCW, 'x-', 'Color',[0.85     0.325	0.098]);
plot(n, Opt_OCW2, 'o-', 'Color',[0        0.447   0.741]);
grid on;
xlabel('number of contending STAs');
ylabel('optimal OCW value');
legend('smallest', 'Largest');
hold off;


figure;
hold on;
plot(n, Opt_EFF1, 'x-', 'Color',[0.85     0.325	0.098]);
plot(n, Opt_EFF2, 'o-', 'Color',[0        0.447   0.741]);
grid on;
xlabel('number of contending STAs');
ylabel('Efficiency');
legend('smallest', 'Largest');
hold off;

Opt_OCW2 = Opt_OCW2';
Opt_EFF2 = Opt_EFF2';

figure;
hold on;
plot(n, Opt_Pc1, 'x-', 'Color',[0.85     0.325	0.098]);
plot(n, Opt_Pc2, 'o-', 'Color',[0        0.447   0.741]);
grid on;
xlabel('number of contending STAs');
ylabel('Collision Rate');
legend('smallest', 'Largest');
hold off;

figure;
hold on;
plot(n, Opt_tau1, 'x-', 'Color',[0.85     0.325	0.098]);
plot(n, Opt_tau2, 'o-', 'Color',[0        0.447   0.741]);
grid on;
xlabel('number of contending STAs');
ylabel('transmission prob.');
legend('smallest', 'Largest');
hold off;

figure;
hold on;
plot(n, Opt_Th, 'x-', 'Color',[0.85     0.325	0.098]);
grid on;
xlabel('number of contending STAs');
ylabel('Throughput');
hold off;


toc;
% t = t';
% p = p';
% eff = eff';
% Th = Th'*10^-6;
% Ds = Ds';
% D = D';
% R = [n' t p eff Th Ds D];
% 
% figure;
% plot(n', eff);
% grid on;
% xlabel('number of stations');
% ylabel('system efficiency');
% 
% figure;
% plot(n', Th);
% grid on;
% xlabel('number of stations');
% ylabel('Throughput(Mb/s)');