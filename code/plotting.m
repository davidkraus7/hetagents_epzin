%% Generate all figures in the paper
% Calls twoassetEZ and twoassetEZ_ep directly.
% Run this script from the code/ directory.

close all
clc

phi = 0;
a_min = 0;
fig_num = 0;

%% Section 5.1: Baseline model (Figures 1-4)

disp('Section 5.1: Baseline model')
[a, c, s, t, gg] = twoassetEZ(4, 0.5, [0.5, 0.5]);

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,c(:,1),'r',a,c(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Consumption','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*2 2])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$c_{L}(a)$', '$c_{H}(a)$'}, ...
       'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,s(:,1),'r',a,s(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Savings','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
yline(0, 'k--', 'LineWidth', 1.5);
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$s_{L}(a)$', '$s_{H}(a)$'}, ...
       'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,t(:,1),'r',a,t(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Risky asset share','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$\theta_{L}(a)$', '$\theta_{H}(a)$'}, ...
       'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,gg(:,1),'r',a,gg(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Density','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$g_{L}(a)$', '$g_{H}(a)$'}, ...
       'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', 30)

%% Section 5.2: Coefficient of relative risk aversion (Figures 5-8)

disp('Section 5.2: Risk aversion comparative statics')
[~, c1, s1, t1, gg1] = twoassetEZ(4, 0.5, [0.5, 0.5]);
[~, c2, s2, t2, gg2] = twoassetEZ(5, 0.5, [0.5, 0.5]);

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,c1(:,1),'r--',a,c1(:,2),'b--',a,c2(:,1),'r',a,c2(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Consumption','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*2 2])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$c_{L,\gamma_-}$', '$c_{H,\gamma_-}$', '$c_{L,\gamma_+}$', '$c_{H,\gamma_+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,s1(:,1),'r--',a,s1(:,2),'b--',a,s2(:,1),'r',a,s2(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Savings','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
yline(0, 'k--', 'LineWidth', 1);
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$s_{L,\gamma_-}$', '$s_{H,\gamma_-}$', '$s_{L,\gamma_+}$', '$s_{H,\gamma_+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,t1(:,1),'r--',a,t1(:,2),'b--',a,t2(:,1),'r',a,t2(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Risky asset share','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$\theta_{L,\gamma_-}$', '$\theta_{H,\gamma_-}$', '$\theta_{L,\gamma_+}$', '$\theta_{H,\gamma_+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,gg1(:,1),'r--',a,gg1(:,2),'b--',a,gg2(:,1),'r',a,gg2(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Density','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$g_{L,\gamma_-}$', '$g_{H,\gamma_-}$', '$g_{L,\gamma_+}$', '$g_{H,\gamma_+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

%% Section 5.3: Intertemporal elasticity of substitution (Figures 9-12)

disp('Section 5.3: IES comparative statics')
[~, c1, s1, t1, gg1] = twoassetEZ(4, 0.4, [0.5, 0.5]);
[~, c2, s2, t2, gg2] = twoassetEZ(4, 0.6, [0.5, 0.5]);

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,c1(:,1),'r--',a,c1(:,2),'b--',a,c2(:,1),'r',a,c2(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Consumption','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*2 2])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$c_{L,\psi_-}$', '$c_{H,\psi_-}$', '$c_{L,\psi_+}$', '$c_{H,\psi_+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,s1(:,1),'r--',a,s1(:,2),'b--',a,s2(:,1),'r',a,s2(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Savings','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
yline(0, 'k--', 'LineWidth', 1);
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$s_{L,\psi_-}$', '$s_{H,\psi_-}$', '$s_{L,\psi_+}$', '$s_{H,\psi_+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,t1(:,1),'r--',a,t1(:,2),'b--',a,t2(:,1),'r',a,t2(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Risky asset share','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$\theta_{L,\psi_-}$', '$\theta_{H,\psi_-}$', '$\theta_{L,\psi_+}$', '$\theta_{H,\psi_+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,gg1(:,1),'r--',a,gg1(:,2),'b--',a,gg2(:,1),'r',a,gg2(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Density','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$g_{L,\psi_-}$', '$g_{H,\psi_-}$', '$g_{L,\psi_+}$', '$g_{H,\psi_+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

%% Section 5.4: Labour income risk (Figures 13-16)

disp('Section 5.4: Labour income risk comparative statics')
[~, c1, s1, t1, gg1] = twoassetEZ(4, 0.5, [0.5, 0.5]);
[~, c2, s2, t2, gg2] = twoassetEZ(4, 0.5, [0.4, 0.6]);

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,c1(:,1),'r--',a,c1(:,2),'b--',a,c2(:,1),'r',a,c2(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Consumption','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*2 2])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$c_{L,-}$', '$c_{H,-}$', '$c_{L,+}$', '$c_{H,+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,s1(:,1),'r--',a,s1(:,2),'b--',a,s2(:,1),'r',a,s2(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Savings','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
yline(0, 'k--', 'LineWidth', 1);
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$s_{L,-}$', '$s_{H,-}$', '$s_{L,+}$', '$s_{H,+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,t1(:,1),'r--',a,t1(:,2),'b--',a,t2(:,1),'r',a,t2(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Risky asset share','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$\theta_{L,-}$', '$\theta_{H,-}$', '$\theta_{L,+}$', '$\theta_{H,+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,gg1(:,1),'r--',a,gg1(:,2),'b--',a,gg2(:,1),'r',a,gg2(:,2),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Density','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$g_{L,-}$', '$g_{H,-}$', '$g_{L,+}$', '$g_{H,+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

%% Section 5.5: Time-varying equity premium (Figures 17-20)

disp('Section 5.5: Time-varying equity premium')
[~, c, s, t, gg] = twoassetEZ_ep([0.06, 0.07], [0.5, 0.5]);
[~, c1, s1, t1, ~] = twoassetEZ(4, 0.5, [0.5, 0.5]);

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,c(:,1),'r--',a,c(:,2),'b--',a,c(:,3),'r',a,c(:,4),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Consumption','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*2 2])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$c_{L,\mu_-}$', '$c_{H,\mu_-}$', '$c_{L,\mu_+}$', '$c_{H,\mu_+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,s(:,1),'r--',a,s(:,2),'b--',a,s(:,3),'r',a,s(:,4),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Savings','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*3 3])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
yline(0, 'k--', 'LineWidth', 1);
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$s_{L,\mu_-}$', '$s_{H,\mu_-}$', '$s_{L,\mu_+}$', '$s_{H,\mu_+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,t(:,1),'r--',a,t(:,2),'b--',a,t(:,3),'r',a,t(:,4),'b', ...
          a,t1(:,1),'r:',a,t1(:,2),'b:','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Risky asset share','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*3 3])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$\theta_{L,\mu_-}$', '$\theta_{H,\mu_-}$', '$\theta_{L,\mu_+}$', '$\theta_{H,\mu_+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

fig_num = fig_num + 1; figure(fig_num)
h1 = plot(a,gg(:,1),'r--',a,gg(:,2),'b--',a,gg(:,3),'r',a,gg(:,4),'b','LineWidth',2);
set(gca,'FontSize',16)
xlabel('Wealth','interpreter','latex', 'FontSize', 35)
ylabel('Density','interpreter','latex', 'FontSize', 35)
xlim([a_min-0.05*5 5])
line([-phi -phi], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
set(gca, 'XTick', []); set(gca, 'YTick', []);
legend(h1, {'$g_{L,\mu_-}$', '$g_{H,\mu_-}$', '$g_{L,\mu_+}$', '$g_{H,\mu_+}$'}, ...
       'Location', 'SouthEast', 'Interpreter', 'latex', 'FontSize', 30)

disp('All figures generated.')
