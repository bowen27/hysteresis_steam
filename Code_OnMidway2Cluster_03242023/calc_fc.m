% Fc = 50:2000:1e6;
% k3 = 0.08;
% fc = k3*log(Fc+1);
% 
% plot(Fc,fc,'-*')

% T = 800:1:1000;
% Tref1 = 200;
% Tref2 = 280;
% a_0 = 0.9;
% k_a = 1e3;
% a = a_0 * exp(-k_a * (1./T - 1/Tref2));
% 
% plot(T, a)

% sigma = 5.67e-8;
% plot(T,sigma.*T.^4)

% hold on
% T = [940 1000 1060];
% E = [300 400 550];
% plot(T,E,'*')
% T1 = 940:1:2000;
% A = 5261.11;
% B = -11.81;
% C = 0.006944;
% E1 = A+B.*T1+C.*T1.^2;
% plot(T1,E1,'--')

T = 200:1:280;
Tref2 = 280;
a_0 = 0.9;
hold on
k_a = 1e3;
a1 = a_0 * exp(-k_a * (1./T - 1/Tref2));
plot(T, a1,'-k','linewidth',2)
k_a = 2e3;
a2 = a_0 * exp(-k_a * (1./T - 1/Tref2));
plot(T, a2,'--k','linewidth',2)
k_a = 3e3;
a3 = a_0 * exp(-k_a * (1./T - 1/Tref2));
plot(T, a3,'-.k','linewidth',2)
k_a = 1e4;
a4 = a_0 * exp(-k_a * (1./T - 1/Tref2));
plot(T, a4,'-ok','linewidth',2)
set(gca,'fontsize',20)
xlabel('temperature (K)')
ylabel('atmos absorptivity')
xlim([200 280])
legend('k_a = 1000','k_a = 2000',...
    'k_a = 3000','k_a = 10000','location','northwest')
