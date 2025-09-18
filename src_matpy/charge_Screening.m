% Charge screening
close all;
clear all;
clc;

clr_arr      = {'g','b','r','m','c',"#ffa500",'k',"#f4a460"};
q1 = 1; eps = [6;25;60;90];
q2 = -0.1;
r = 2:0.1:10;
figure
hold on
box on
for i = 1:length(eps)
  uofr = (q1*q2)./(eps(i).*r);
  plot(r, uofr, 'LineWidth', 2, 'color',clr_arr{i},'LineWidth',2);
  leg_arr{i} = ['\epsilon:' num2str(eps(i))];
end
xlabel('r','FontSize',20);
ylabel('u(r)','FontSize',20);
legend(leg_arr,'location','southeast','FontSize',20)
set(gca, 'FontSize', 16)

leg_arr = {}
eps_ref = 1;
q2 = -0.05:-0.2:-1;
figure
hold on
box on
for i = 1:length(q2)
  uofr = (q1.*q2(i))./(eps_ref.*r);
  plot(r, uofr, 'LineWidth', 2, 'color',clr_arr{i},'LineWidth',2);
  leg_arr{i} = ['q: ' num2str(q2(i))];
end
xlabel('r','FontSize',20);
ylabel('u(r)','FontSize',20);
legend(leg_arr,'location','southeast','FontSize',20)
set(gca, 'FontSize', 16)

