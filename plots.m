set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% Van der Pol
figure(1000)
for i = 1:length(S_van)
hold on
plot(Polyhedron(S_van{i}))
end

xlabel("$x_1$")
ylabel("$x_2$")
title("Van der Pol Oscillator")
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4 3]);
set(gcf, 'PaperSize', [4 3]);
saveas(gcf,'img/van_der_pol.pdf')

%% cartpole pendulum
figure(1001)
for i = 1:length(S_cartpole_pendulum)
hold on
plot(Polyhedron(S_cartpole_pendulum{i}))
end

xlabel("$x_1$")
ylabel("$x_2$")
title("Cartpole (pendulum only)")
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4 3]);
set(gcf, 'PaperSize', [4 3]);
saveas(gcf,'img/cartpole_pendulum.pdf')