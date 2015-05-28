function plotBodyCirculation (inputFile)

system('./main 1');

close all

X = importdata(inputFile);

t = X(:,1);
dt = t(2) - t(1);

Qao = X(:,12);
Qtr = X(:,20);
Qpv = X(:,21);
Qmi = X(:,29);

Vlv = 111.8*ones(length(t));
Vrv = 108.6*ones(length(t));

for i = 2:length(t);
    Vlv(i) = Vlv(i-1) + (Qmi(i) - Qao(i) + X(i,31)) * dt;
    Vrv(i) = Vrv(i-1) + (Qtr(i) - Qpv(i) + X(i,32)) * dt;
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
for i = 12:2:20 plot(X(:,1), X(:,i)); end
plot(t, X(:,30))
legend('Q_{ao}', 'Q_{sa}', 'Q_{sp}', 'Q_{sv}', 'Q_{tr}')
xlabel('t [s]')
ylabel('Q [ml/s]')
%axis([t1 t2 -100 250])
%xlim([t1 t2])
grid on

subplot(2,3,3)
hold on
for i = 32:8:40 plot(X(:,1), X(:,i)); end
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
for i = 21:2:29 plot(X(:,1), X(:,i)); end
plot(t, X(:,31))
legend('Q_{pv}', 'Q_{pa}', 'Q_{pp}', 'Q_{pv}', 'Q_{mi}')
xlabel('t [s]')
ylabel('Q [ml/s]')
%axis([t1 t2 -100 250])
%xlim([t1 t2])
grid on

subplot(2,3,6)
hold on
for i = 41:8:49 plot(X(:,1), X(:,i)); end
legend('N_{pv}', 'N_{mi}')
xlabel('t [s]')
ylabel('N')
%axis([t1 t2 0 1])
grid on


end