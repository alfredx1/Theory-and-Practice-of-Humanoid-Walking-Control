
load alfred_hw6.txt;

x=alfred_hw6(:,1); 
y1=alfred_hw6(:,2);
y2=alfred_hw6(:,3);
y3=alfred_hw6(:,4);
y4=alfred_hw6(:,5);
y5=alfred_hw6(:,6);

y6=alfred_hw6(:,7);


%% X
figure (1)
plot(x,y1);
xlim([0, 11]);
hold on;
plot(x,y3);
hold on;
plot(x,y5, '--');
ylim([-0.1, 0.1]);
hold on;
xlabel('Time [sec]');
ylabel('Displacement [m]');
title('X motion');
legend('ZMP','CoM', 'Calculated ZMP')
hold on;
grid on;

%%
figure (2)
plot(x,y2);
xlim([0, 11]);
hold on;
plot(x,y4);
hold on;
plot(x,y6,'--');
ylim([-0.15, 0.15]);
xlabel('Time [sec]');
ylabel('Displacement [m]');
title('Y motion');
legend('ZMP','CoM', 'Calculated ZMP')
hold on;
grid on;

