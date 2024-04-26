figure(2);
subplot(3,2,1);
hold on;
plot((NTS-NTSPP:NTS)*dt/(1/freq),CD(NTS-NTSPP:NTS),'black','LineWidth',1);
xlabel('$t/T$','interpreter','latex','fontsize',15);
ylabel('$C_d$','interpreter','latex','fontsize',15,'rotation',0);

subplot(3,2,2);
hold on;
plot(ALPHA_EFF(NTS-NTSPP:NTS),CD(NTS-NTSPP:NTS),'*');
xlabel('$\alpha_{eff}$','interpreter','latex','fontsize',15);
ylabel('$C_d$','interpreter','latex','fontsize',15,'rotation',0);

subplot(3,2,3);
hold on;
plot((NTS-NTSPP:NTS)*dt/(1/freq),CL(NTS-NTSPP:NTS),'black','LineWidth',1);
xlabel('$t/T$','interpreter','latex','fontsize',15);
ylabel('$C_l$','interpreter','latex','fontsize',15,'rotation',0);

subplot(3,2,4);
hold on;
plot(ALPHA_EFF(NTS-NTSPP:NTS),CL(NTS-NTSPP:NTS),'*');
xlabel('$\alpha_{eff}$','interpreter','latex','fontsize',15);
ylabel('$C_l$','interpreter','latex','fontsize',15,'rotation',0);

subplot(3,2,5);
hold on;
plot((NTS-NTSPP:NTS)*dt/(1/freq),CM(NTS-NTSPP:NTS),'black','LineWidth',1);
xlabel('$t/T$','interpreter','latex','fontsize',15);
ylabel('$C_m$','interpreter','latex','fontsize',15,'rotation',0);

subplot(3,2,6);
hold on;
plot(ALPHA_EFF(NTS-NTSPP:NTS),CM(NTS-NTSPP:NTS),'*');
xlabel('$\alpha_{eff}$','interpreter','latex','fontsize',15);
ylabel('$C_m$','interpreter','latex','fontsize',15,'rotation',0);

pause(0.01);