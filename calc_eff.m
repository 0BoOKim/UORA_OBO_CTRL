clear all;
close all;

n = 1:200;%[1 5:5:200];      % number of STAs
W0 = 31;         % OCWmin
m = 5;          % number of backoff stages
M = 8;          % number of RA-RUs 
L_MPDU = 2000;  % Byte
Data_Rate = 6.6667*10^6;%65*10^6; %;  

L_PHY = 40*10^-6;      % preamble length (40 usec)
L_Trigger = 100*10^-6; % length of trigger frame
L_BACK = 68*10^-6;     % length of block ack
SIFS = 16*10^-6;

T = 8*L_MPDU/Data_Rate + 3*SIFS + 3*L_PHY + L_BACK;


for i = 1:length(n)
    t(i)=fzero(@tau,[0,1],[],n(i),W0,m,M);   % tau
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
t = t';
p = p';
eff = eff';
Th = Th'*10^-6;
Ds = Ds';
D = D';
R = [n' t p eff Th Ds D];

figure;
plot(n', eff);
grid on;
xlabel('number of stations');
ylabel('system efficiency');

figure;
plot(n', Th);
grid on;
xlabel('number of stations');
ylabel('Throughput(Mb/s)');