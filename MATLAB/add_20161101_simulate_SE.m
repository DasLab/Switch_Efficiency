R = 8.314;
RNA = [1:100] * 1e-9;
T = [0:40] + 273.15;
L_high = 1.5 * 1e-6;
L_low = 0.01 * 1e-6;

%%% 2000 equiv Mg2+
dHd = 110 * 1e3;
dSd = 238;
% dGd = 39 * 1e3;
%%% SI fig 7
dHpre = 59 * 1e3;
dSpre = 187;
% dGpre = 2 * 1e3;

% grid search for RNA conc. vs. T
SE_3 = zeros(length(RNA), length(T));
SE_2 = zeros(length(RNA), length(T));
SE_3as2 = zeros(length(RNA), length(T));
for j = 1:length(T);
    dGd = dHd - dSd * T(j);
    dGpre = dHpre - dSpre * T(j);
    Kd = exp(-dGd / (R * T(j)));
    Kpre = exp(-dGpre / (R * T(j)));
    for i = 1:length(RNA);
        [holo_high, apoA_high, ~] = add_fractions_3state(Kd, Kpre, L_high, RNA(i));
        [holo_low, apoA_low, ~] = add_fractions_3state(Kd, Kpre, L_low, RNA(i));
        SE_3(i,j) = (holo_high - holo_low) / RNA(i);
        SE_3as2(i,j) = SE_3(i,j) + (apoA_high - apoA_low) / RNA(i);
        
        [holo_high, ~] = add_fractions_2state(Kd, L_high, RNA(i));
        [holo_low, ~] = add_fractions_2state(Kd, L_low, RNA(i));
        SE_2(i,j) = (holo_high - holo_low) / RNA(i);
    end;
end;

figure();
subplot(1,3,1);
[c, h] = contourf(SE_3, 0:0.1:1); clabel(c), caxis([0 1]); colormap(parula);
title('3-state (apoA = OFF)', 'fontsize',20,'fontweight','bold');
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[RNA] (nM)', 'fontsize',16,'fontweight','bold');
set(gca,'ydir','normal', 'xtick',[1:10:length(T)], 'xticklabel', 0:10:length(T)); colorbar;
subplot(1,3,2);
[c, h] = contourf(SE_2, 0:0.1:1); clabel(c), caxis([0 1]); colormap(parula);
title('2-state', 'fontsize',20,'fontweight','bold');
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[RNA] (nM)', 'fontsize',16,'fontweight','bold');
set(gca,'ydir','normal', 'xtick',[1:10:length(T)], 'xticklabel', 0:10:length(T)); colorbar;
subplot(1,3,3);
[c, h] = contourf(SE_3as2, 0:0.1:1); clabel(c), caxis([0 1]); colormap(parula);
title('3-state (apoA = ON)', 'fontsize',20,'fontweight','bold');
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[RNA] (nM)', 'fontsize',16,'fontweight','bold');
set(gca,'ydir','normal', 'xtick',[1:10:length(T)], 'xticklabel', 0:10:length(T)); colorbar;


%%% single plots
% use RNA = 20 nM
holo_high_20 = zeros(1, length(T));
holo_low_20 = zeros(1, length(T));
apoA_high_20 = zeros(1, length(T));
apoA_low_20 = zeros(1, length(T));
holo_high_std = zeros(1, length(T));
holo_low_std = zeros(1, length(T));
RNA_20 = 20 * 1e-9;
L_high = 1.0 * 1e-6;

for j = 1:length(T);
    dGd = dHd - dSd * T(j);
    dGpre = dHpre - dSpre * T(j);
    Kd = exp(-dGd / (R * T(j)));
    Kpre = exp(-dGpre / (R * T(j)));
    [holo_high_20(j), apoA_high_20(j), ~] = add_fractions_3state(Kd, Kpre, L_high, RNA_20);
    [holo_low_20(j), apoA_low_20(j), ~] = add_fractions_3state(Kd, Kpre, L_low, RNA_20);
    [holo_high_std(j), ~] = add_fractions_2state(Kd, L_high, RNA_20);
    [holo_low_std(j), ~] = add_fractions_2state(Kd, L_low, RNA_20);
end;

holo_high_20 = holo_high_20 / RNA_20;
holo_low_20 = holo_low_20 / RNA_20;
apoA_high_20 = apoA_high_20 / RNA_20;
apoA_low_20 = apoA_low_20 / RNA_20;
holo_high_std = holo_high_std / RNA_20;
holo_low_std = holo_low_std / RNA_20;

figure();
subplot(3,1,1);
plot(holo_high_20, 'bo-', 'linewidth',2); hold on;
plot(holo_low_20, 'co-', 'linewidth',2);
title('3-state (apoA = OFF) SE', 'fontsize',20,'fontweight','bold'); axis([1 45 0 1]);
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[holo] (%)', 'fontsize',16,'fontweight','bold');
set(gca,'xtick', 0:5:length(T), 'xticklabel', 0:5:length(T), 'fontsize',16); legend('1.0 uM adenine', '0.01 uM adenine', 'location', 'northeast');
subplot(3,1,2);
plot(holo_high_std, 'ro-', 'linewidth',2); hold on;
plot(holo_low_std, 'mo-', 'linewidth',2);
title('2-state SE', 'fontsize',20,'fontweight','bold'); axis([1 45 0 1]);
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[holo] (%)', 'fontsize',16,'fontweight','bold');
set(gca,'xtick', 0:5:length(T),'xticklabel', 0:5:length(T), 'fontsize',16); legend('1.0 uM adenine', '0.01 uM adenine', 'location', 'northeast');
subplot(3,1,3);
plot(holo_high_20 + apoA_high_20, 'ko-', 'linewidth',2); hold on;
plot(holo_low_20 + apoA_low_20, 'o-', 'color', [0.5 0.5 0.5], 'linewidth',2);
title('3-state (apoA = ON) SE', 'fontsize',20,'fontweight','bold'); axis([1 45 0 1]);
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[holo] + [apoA] (%)', 'fontsize',16,'fontweight','bold');
set(gca,'xtick', 0:5:length(T),'xticklabel', 0:5:length(T), 'fontsize',16); legend('1.0 uM adenine', '0.01 uM adenine', 'location', 'northeast');

%%%
figure();
subplot(3,2,1);
plot(holo_high_20, 'bo-', 'linewidth',2); hold on;
plot(holo_low_20, 'co-', 'linewidth',2);
make_lines(29.5,'g',2);
title('3-state (apoA = OFF) SE', 'fontsize',20,'fontweight','bold'); axis([1 40 0 1]);
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[holo] (%)', 'fontsize',16,'fontweight','bold');
set(gca,'xtick', 0:5:length(T), 'xticklabel', 0:5:length(T), 'fontsize',16); legend('1.0 uM adenine', '0.01 uM adenine', 'location', 'northeast');
subplot(3,2,2);
plot(holo_high_std, 'ro-', 'linewidth',2); hold on;
plot(holo_low_std, 'mo-', 'linewidth',2);
make_lines(29.5,'g',2);
title('2-state SE', 'fontsize',20,'fontweight','bold'); axis([1 40 0 1]);
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[holo] (%)', 'fontsize',16,'fontweight','bold');
set(gca,'xtick', 0:5:length(T),'xticklabel', 0:5:length(T), 'fontsize',16); legend('1.0 uM adenine', '0.01 uM adenine', 'location', 'east');
subplot(3,2,[3 5]);
[c, h] = contourf(SE_3, 0:0.1:1); clabel(c), caxis([0 1]); colormap(parula);
title('3-state (apoA = OFF)', 'fontsize',20,'fontweight','bold');
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[RNA] (nM)', 'fontsize',16,'fontweight','bold');
set(gca,'ydir','normal', 'xtick',[1:10:length(T)], 'xticklabel', 0:10:length(T));
subplot(3,2,[4 6]);
[c, h] = contourf(SE_2, 0:0.1:1); clabel(c), caxis([0 1]); colormap(parula);
title('2-state', 'fontsize',20,'fontweight','bold');
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[RNA] (nM)', 'fontsize',16,'fontweight','bold');
set(gca,'ydir','normal', 'xtick',[1:10:length(T)], 'xticklabel', 0:10:length(T)); 

figure();
subplot(3,2,1);
plot(holo_high_20 + apoA_high_20, 'ko-', 'linewidth',2); hold on;
plot(holo_low_20 + apoA_low_20, 'o-', 'color', [0.5 0.5 0.5], 'linewidth',2);
make_lines(29.5,'g',2);
title('3-state (apoA = ON) SE', 'fontsize',20,'fontweight','bold'); axis([1 40 0 1]);
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[holo] + [apoA] (%)', 'fontsize',16,'fontweight','bold');
set(gca,'xtick', 0:5:length(T),'xticklabel', 0:5:length(T), 'fontsize',16); legend('1.0 uM adenine', '0.01 uM adenine', 'location', 'northeast');
subplot(3,2,[3 5]);
[c, h] = contourf(SE_3as2, 0:0.1:1); clabel(c), caxis([0 1]); colormap(parula);
title('3-state (apoA = ON)', 'fontsize',20,'fontweight','bold');
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[RNA] (nM)', 'fontsize',16,'fontweight','bold');
set(gca,'ydir','normal', 'xtick',[1:10:length(T)], 'xticklabel', 0:10:length(T));
