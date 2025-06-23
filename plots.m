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

x_0 = [-2;0];
robotSize = 0.25
X_0 = interval(-robotSize*[0.75;1], robotSize/2*[0.75;1]).mptPolytope.P + x_0

vidfile = VideoWriter('robot.avi')
open(vidfile)
for t = -pi+0.001:pi/64:pi
    t
    R = [cos(t), sin(t); -sin(t), cos(t)]
    X  = R*X_0
    clf; hold on;
    for i = 1:length(S_p)
        tmp = S_p{i}(3,:);
        if tmp.contains(t)
            plot(S_p{i}.mptPolytope.P.projection(1:2), 'color', 'r')
        end
    end
    plot(X, 'color', 'k');
    axis equal
    getframe(gcf)
    cdata = print('-RGBImage','-r200','-noui');
    frame = im2frame(cdata);
    writeVideo(vidfile,frame)
end
close(vidfile)

%%
figure
S_p = importfile3d("data_exploration/s13.txt");
tiledlayout(1,2,'TileSpacing','compact', 'Padding','compact')

k = 1;
angles = {'-\frac{3\pi}{4}', '-\frac{\pi}{4}'};
for t = [-3*pi/4, -pi/4]
    nexttile
    hold on;
    for i = 1:length(S_p)
        tmp = S_p{i}(3,:);
        if tmp.contains(t)
            plot(S_p{i}.mptPolytope.P.projection(1:2), 'color', 'r')
        end
    end
    axis equal
    title(strcat("$\theta = ", angles{k}, "$"))
    xlabel('$x_1$')
    ylabel('$x_2$')
    k = k + 1;
end

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6 3]);
set(gcf, 'PaperSize', [6 3]);
saveas(gcf,'img/exploration_slice.pdf')
%%
figure
S_p = importfile3d("data_exploration/s13.txt");
hold on
for i = 1:length(S_p)
    plot(S_p{i}.mptPolytope.P, 'color', 'r')
end

xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$\theta$')

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4 3]);
set(gcf, 'PaperSize', [4 3]);
saveas(gcf,'img/exploration_3d.pdf')

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