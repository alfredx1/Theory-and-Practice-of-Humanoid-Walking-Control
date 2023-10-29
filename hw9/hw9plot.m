
load alfred_hw9.txt;
load alfred_hw9_nocpt.txt;

x=alfred_hw9(:,1);
y1=alfred_hw9(:,2);
y2=alfred_hw9(:,3);
y3=alfred_hw9(:,4);
y4=alfred_hw9(:,5);

x1=alfred_hw9_nocpt(:,1);
y5=alfred_hw9_nocpt(:,2);
y6=alfred_hw9_nocpt(:,3);
y7=alfred_hw9_nocpt(:,4);
y8=alfred_hw9_nocpt(:,5);

%% X %%
figure (1)
plot(x,y1);
xlim([0, 11]);
hold on;
plot(x,y3);
hold on;
plot(x1,y5);
hold on;
plot(x1,y7);
xlabel('Time [sec]');
ylabel('Displacement [m]');
title('X motion');
legend('CP_px (CPT)', 'cpx(CPT)', 'CP_px', 'cpx(CPT)');
hold on;
grid on;

%% Y %%
figure (2)
plot(x,y2);
xlim([0, 11]);
hold on;
plot(x,y4);
hold on;
plot(x1,y6);
hold on;
plot(x1,y8);
hold on;
xlabel('Time [sec]');

ylabel('Displacement [m]');
title('Y motion');
legend('CP_py (CPT)', 'cpy(CPT)', 'CP_py', 'cpy(CPT)');
hold on;
grid on;

