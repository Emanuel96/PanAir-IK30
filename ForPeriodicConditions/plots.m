figure(2);
subplot(3,2,1);
plot((1:k)*dt/(1/freq),CD(1:k),'black','LineWidth',1);
xlabel('$t/T$','interpreter','latex','fontsize',15);
ylabel('$C_d$','interpreter','latex','fontsize',15,'rotation',0);

subplot(3,2,2);
hold on;
plot(ALPHA_EFF(1:k),CD(1:k),'*');
xlabel('$\alpha_{eff}$','interpreter','latex','fontsize',15);
ylabel('$C_d$','interpreter','latex','fontsize',15,'rotation',0);

subplot(3,2,3);
plot((1:k)*dt/(1/freq),CL(1:k),'black','LineWidth',1);
xlabel('$t/T$','interpreter','latex','fontsize',15);
ylabel('$C_l$','interpreter','latex','fontsize',15,'rotation',0);

subplot(3,2,4);
hold on;
plot(ALPHA_EFF(1:k),CL(1:k),'*');
xlabel('$\alpha_{eff}$','interpreter','latex','fontsize',15);
ylabel('$C_l$','interpreter','latex','fontsize',15,'rotation',0);

subplot(3,2,5);
plot((1:k)*dt/(1/freq),CM(1:k),'black','LineWidth',1);
xlabel('$t/T$','interpreter','latex','fontsize',15);
ylabel('$C_m$','interpreter','latex','fontsize',15,'rotation',0);

subplot(3,2,6);
hold on;
plot(ALPHA_EFF(1:k),CM(1:k),'*');
xlabel('$\alpha_{eff}$','interpreter','latex','fontsize',15);
ylabel('$C_m$','interpreter','latex','fontsize',15,'rotation',0);