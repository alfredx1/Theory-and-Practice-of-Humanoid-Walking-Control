
load alfred_hw7.txt;

x=alfred_hw7(:,1);
y1=alfred_hw7(:,2);
y2=alfred_hw7(:,3);
y3=alfred_hw7(:,4);
y4=alfred_hw7(:,5);
y5=alfred_hw7(:,6);
y6=alfred_hw7(:,7);


%% X
figure (1)
plot(x,y5);
xlim([0, 11]);
hold on;
plot(x,y1);
hold on;
plot(x,y3, '--');
ylim([-0.1, 0.1]);
hold on;
xlabel('Time [sec]');
ylabel('Displacement [m]');
title('X motion');
legend('Reference ZMP','Reference COM', 'Online COM');
hold on;
grid on;

%%
figure (2)
plot(x,y6);
xlim([0, 11]);
hold on;
plot(x,y2);
hold on;
plot(x,y4,'--');
ylim([-0.15, 0.15]);
xlabel('Time [sec]');
ylabel('Displacement [m]');
title('Y motion');
legend('Reference ZMP','Reference COM', 'Online COM');
hold on;
grid on;

