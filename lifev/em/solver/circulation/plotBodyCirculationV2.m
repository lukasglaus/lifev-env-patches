function plotBodyCirculationV2 (inputFile)

system('./main 1 0.0025');

close all

X = importdata(inputFile);

t = X(:,1);
dt = t(2) - t(1);

Qao = X(:,12);
Qtr = X(:,15);
Qpv = X(:,13);
Qmi = X(:,14);

Vlv = 111.8*ones(length(t));
Vrv = 108.6*ones(length(t));

for i = 2:length(t);
    Vlv(i) = Vlv(i-1) + (Qmi(i) - Qao(i)) * dt;
    Vrv(i) = Vrv(i-1) + (Qtr(i) - Qpv(i)) * dt;
end

figure(2)
plot( t, Vlv, t, Vrv )
legend('LV', 'RV')
xlabel('t [s]')
ylabel('Volume [ml]')
grid on

figure(1)
subplot(2,3,1)
hold on
for i = 2:6 plot(X(:,1), X(:,i)); end
legend('p_{lv}','p_{sa}', 'p_{sp}', 'p_{sv}', 'p_{ra}')
xlabel('t [s]')
ylabel('p [mmHg]')
%axis([t1 t2 0 16])
%xlim([t1 t2])
grid on

subplot(2,3,2)
hold on
for i = 12:2:16 plot(X(:,1), X(:,i)); end
%plot(t, X(:,30))
legend('Q_{ao}', 'Q_{sa}', 'Q_{sp}', 'Q_{sv}', 'Q_{tr}')
xlabel('t [s]')
ylabel('Q [ml/s]')
%axis([t1 t2 -100 250])
%xlim([t1 t2])
grid on


Nao = X(:,12)*3.75e-4.*X(:,2)./( X(:,2) - X(:,3) );
Ntr = X(:,15)*2.5e-3./( X(:,11) - X(:,2) );

subplot(2,3,3)
hold on
plot(t, Nao, t, Ntr)
legend('N_{ao}', 'N_{tr}')
xlabel('t [s]')
ylabel('N')
%axis([t1 t2 0 1])
grid on

subplot(2,3,4)
hold on
for i = 7:11 plot(X(:,1), X(:,i)); end
legend('p_{rv}','p_{pa}', 'p_{pp}', 'p_{pv}', 'p_{la}')
xlabel('t [s]')
ylabel('p [mmHg]')
%axis([t1 t2 0 16])
%xlim([t1 t2])
grid on

subplot(2,3,5)
hold on
for i = 13:2:17 plot(X(:,1), X(:,i)); end
%plot(t, X(:,31))
legend('Q_{pv}', 'Q_{pa}', 'Q_{pp}', 'Q_{pv}', 'Q_{mi}')
xlabel('t [s]')
ylabel('Q [ml/s]')
%axis([t1 t2 -100 250])
%xlim([t1 t2])
grid on

Npv = X(:,13)*1.4e-3.*X(:,7)./( X(:,7) - X(:,8) );
Nmi = X(:,14)*2.5e-3./( X(:,6) - X(:,7) );

subplot(2,3,6)
hold on
plot(t, Npv, t, Nmi)
legend('N_{pv}', 'N_{mi}')
xlabel('t [s]')
ylabel('N')
%axis([t1 t2 0 1])
grid on


end