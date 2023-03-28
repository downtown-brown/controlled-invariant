set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  

figure(1)
hold on
plot(Polyhedron(P{1}))
plot(S{i}.mptPolytope.P, 'color', 'blue')
for i = 1:size(P,1)
plot(Polyhedron(P{i}));
end

for i = 1:size(S,1)
plot(S{i}.mptPolytope.P, 'color', 'blue');
end

l = legend('Ours', 'x');
l.Location = 'southeast';
xlabel('$x_1$')
ylabel('$x_2$')

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 5 4]);
saveas(gcf,'img/x.eps','epsc')
saveas(gcf,'img/x.png')

%%
figure(2)
U2 = {};
for i = 1:size(U,1)
   U2{i} = vertcat(U{i}{:});
end
for i = 1:size(P,1)
    tmp = interval([min(P{i}(:,1)); min(P{i}(:,2)); min(U2{i}(:,2))/0.5], [max(P{i}(:,1)); max(P{i}(:,2)); max(U2{i}(:,2))/0.5]);
plot(tmp.mptPolytope.P);
hold on
end

xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$\mathcal{U}([x])$')

%%
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 5 3]);
saveas(gcf,'img/e.eps','epsc')
saveas(gcf,'img/e.png')
