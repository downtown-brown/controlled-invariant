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

%%
figure
S_p = importfile3d("data_exploration/s13.txt");

for t = -pi/2:pi/4:pi/2
    clf; hold on;
    for i = 1:length(S_p)
        tmp = S_p{i}(3,:);
        if tmp.contains(t)
            plot(S_p{i}.mptPolytope.P, 'color', 'r')
        end
    end
    view(0,90)
    axis equal
end

%%
figure
S_p = importfile3d("data_exploration/s12.txt");
for i = 1:length(S_p)
    tmp = S_p{i}(3,:);
    plot(S_p{i}.mptPolytope.P, 'color', 'r')
    hold on
end

%%
S_p = importfile2d("data_cart_ch/s3.txt");
for i = 1:length(S_p)
hold on
tmp = S_p{i};
plot(S_p{i}.mptPolytope.P, 'color', 'r', 'alpha', 0.5)
end

% xlabel("$x_1$")
% ylabel("$x_2$")
% title("Cartpole (pendulum only)")
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0 0 4 3]);
% set(gcf, 'PaperSize', [4 3]);
% saveas(gcf,'img/cartpole_pendulum.pdf')

%%

figure
for i = 1:length(S)
hold on
tmp = S{i};
tmpi = interval(tmp(:,1), tmp(:,2));
plot(tmpi.mptPolytope.P.projection([2,3,4]))
end

% xlabel("$x_2$")
% ylabel("$x_3$")
% zlabel("$x_4$")
% title("Cartpole (3D projection)")
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0 0 4 3]);
% set(gcf, 'PaperSize', [4 3]);
% saveas(gcf,'img/cartpole.pdf')