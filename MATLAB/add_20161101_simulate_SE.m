R = 8.31446;
RNA = [1:100] * 1e-9;
T = [0:30] + 273.15;
L_high = 1.5 * 1e-6;
L_low = 0.01 * 1e-6;

%%% 2000 equiv Mg2+
dHd = 110.470 * 1e3;
dSd = 238.303;
% dGd = 39 * 1e3;
%%% SI fig 7
%dHpre = 59 * 1e3;
%dSpre = 187;
dHpre = 52.0028 * 1e3;  % from Furtig .nb
dSpre = 167.281; % from Furtig .nb
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

set( figure(1), 'Position', [37 522 1140 283] );
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
holo_high_test  = zeros(1, length(T));
holo_low_test   = zeros(1, length(T));
apoA_high_test  = zeros(1, length(T));
apoA_low_test   = zeros(1, length(T));
holo_high_std = zeros(1, length(T));
holo_low_std  = zeros(1, length(T));
%RNA_test = 20 * 1e-9;
RNA_test = 1.5 * 1e-9; % from Furtig .nb
L_high = 1.0 * 1e-6;

for j = 1:length(T);
    dGd = dHd - dSd * T(j);
    dGpre = dHpre - dSpre * T(j);
    Kd = exp(-dGd / (R * T(j)));
    Kpre = exp(-dGpre / (R * T(j)));
    [holo_high_test(j), apoA_high_test(j), ~] = add_fractions_3state(Kd, Kpre, L_high, RNA_test);
    [holo_low_test(j),  apoA_low_test(j),  ~] = add_fractions_3state(Kd, Kpre, L_low, RNA_test);
    [holo_high_std(j), ~] = add_fractions_2state(Kd, L_high, RNA_test);
    [holo_low_std(j), ~]  = add_fractions_2state(Kd, L_low, RNA_test);
end;

holo_high_test = holo_high_test / RNA_test;
holo_low_test = holo_low_test / RNA_test;
apoA_high_test = apoA_high_test / RNA_test;
apoA_low_test = apoA_low_test / RNA_test;
holo_high_std = holo_high_std / RNA_test;
holo_low_std = holo_low_std / RNA_test;

set( figure(2), 'Position', [271 212 348 593] );
subplot(3,1,1);
plot(holo_high_test, 'bo-', 'linewidth',2); hold on;
plot(holo_low_test, 'co-', 'linewidth',2);
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
plot(holo_high_test + apoA_high_test, 'ko-', 'linewidth',2); hold on;
plot(holo_low_test + apoA_low_test, 'o-', 'color', [0.5 0.5 0.5], 'linewidth',2);
title('3-state (apoA = ON) SE', 'fontsize',20,'fontweight','bold'); axis([1 45 0 1]);
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[holo] + [apoA] (%)', 'fontsize',16,'fontweight','bold');
set(gca,'xtick', 0:5:length(T),'xticklabel', 0:5:length(T), 'fontsize',16); legend('1.0 uM adenine', '0.01 uM adenine', 'location', 'northeast');

%%%
set( figure(3), 'Position',  [ 44   292   699   482 ] );
subplot(3,2,1);
plot(T, holo_high_test, 'k-', 'linewidth',2); hold on;
plot(T, holo_low_test, 'b-', 'linewidth',2);
title('3-state (apoA = OFF) SE', 'fontsize',20,'fontweight','bold');  axis([273 273+30 0 1]);
xlabel('T (K)', 'fontsize',16,'fontweight','bold'); ylabel('[holo] (%)', 'fontsize',16,'fontweight','bold');
set(gca,'xtick', 275+[0:5:length(T)],'xticklabel', 275+[0:5:length(T)], 'fontsize',16); legend('1.0 uM adenine', '0.01 uM adenine', 'location', 'northeast');

subplot(3,2,2);
plot(T,holo_low_std,'k','linew',2);hold on
plot(T,holo_high_std,'b','linew',2);
title('2-state SE', 'fontsize',20,'fontweight','bold'); axis([273 273+30 0 1]);
xlabel('T (K)', 'fontsize',16,'fontweight','bold'); ylabel('[holo] (%)', 'fontsize',16,'fontweight','bold');
set(gca,'xtick', 275+[0:5:length(T)],'xticklabel', 275+[0:5:length(T)], 'fontsize',16); legend('1.0 uM adenine', '0.01 uM adenine', 'location', 'east');

subplot(3,2,[3 5]);
[c, h] = contourf(T, RNA/1e-9, SE_3, 0:0.1:1); clabel(c), caxis([0 1]); colormap(parula);
title('3-state (apoA = OFF)', 'fontsize',20,'fontweight','bold');
xlabel('T (K)', 'fontsize',16,'fontweight','bold'); ylabel('[RNA] (nM)', 'fontsize',16,'fontweight','bold');
set(gca,'ydir','normal', 'xtick',[1:10:length(T)], 'xticklabel', 0:10:length(T));

subplot(3,2,[4 6]);
[c, h] = contourf(T,RNA/1e-9,SE_2, 0:0.1:1); clabel(c), caxis([0 1]); colormap(parula);
title('2-state', 'fontsize',20,'fontweight','bold');
xlabel('T (K)', 'fontsize',16,'fontweight','bold'); ylabel('[RNA] (nM)', 'fontsize',16,'fontweight','bold');
set(gca,'ydir','normal', 'xtick',275+[0:5:300], 'xticklabel', 275+[0:5:300]); 

%%%
figure(4);
subplot(3,2,1);
plot(holo_high_test + apoA_high_test, 'ko-', 'linewidth',2); hold on;
plot(holo_low_test + apoA_low_test, 'o-', 'color', [0.5 0.5 0.5], 'linewidth',2);
make_lines(29.5,'g',2);
title('3-state (apoA = ON) SE', 'fontsize',20,'fontweight','bold'); axis([1 40 0 1]);
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[holo] + [apoA] (%)', 'fontsize',16,'fontweight','bold');
set(gca,'xtick', 0:5:length(T),'xticklabel', 0:5:length(T), 'fontsize',16); legend('1.0 uM adenine', '0.01 uM adenine', 'location', 'northeast');

subplot(3,2,[3 5]);
[c, h] = contourf(SE_3as2, 0:0.1:1); clabel(c), caxis([0 1]); colormap(parula);
title('3-state (apoA = ON)', 'fontsize',20,'fontweight','bold');
xlabel('T (C)', 'fontsize',16,'fontweight','bold'); ylabel('[RNA] (nM)', 'fontsize',16,'fontweight','bold');
set(gca,'ydir','normal', 'xtick',[1:10:length(T)], 'xticklabel', 0:10:length(T));


figure(3); % the most important plot!