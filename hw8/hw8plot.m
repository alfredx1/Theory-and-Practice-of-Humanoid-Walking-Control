
load alfred_hw8.txt;

x=alfred_hw8(:,1);
y1=alfred_hw8(:,2);
y2=alfred_hw8(:,3);
y3=alfred_hw8(:,4);
y4=alfred_hw8(:,5);
y5=alfred_hw8(:,6);
y6=alfred_hw8(:,7);
y7=alfred_hw8(:,8);
y8=alfred_hw8(:,9);
y9=alfred_hw8(:,10);

%% X %%
figure (1)
plot(x,y7);
xlim([0, 11]);
hold on;
plot(x,y5,'--');
hold on;
plot(x,y3);
ylim([-0.1, 0.1]);
hold on;
plot(x,y1);
hold on;
xlabel('Time [sec]');
ylabel('Displacement [m]');
title('X motion');
legend('CP eos', 'CP', 'ZMP', 'CoM');
hold on;
grid on;

%% Y %%
figure (2)
plot(x,y8);
xlim([0, 11]);
hold on;
plot(x,y6,'--');
hold on;
plot(x,y4);
ylim([-0.15, 0.15]);
hold on;
plot(x,y2);
hold on;
xlabel('Time [sec]');
ylabel('Displacement [m]');
title('Y motion');
legend('CP eos', 'CP', 'ZMP', 'CoM');
hold on;
grid on;


figure(3);
plot(x,y9);
ylabel('Td [sec]');
xlabel('Time [sec}')
