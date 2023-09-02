set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  

figure(1)
hold on
plot(Polyhedron(P{1}), 'color', 'cyan')
plot(S{1}.mptPolytope.P, 'color', [0 0.4470 0.7410])
for i = 1:size(P,1)
plot(Polyhedron(P{i}), 'color', 'cyan');
end

for i = 1:size(S,1)
plot(S{i}.mptPolytope.P, 'color', [0 0.4470 0.7410]);
end

l = legend('Algs 2, 3, \& [9] $n_u=1000$', '[9] $n_u=10$');
l.Location = 'northoutside';
xlabel('$x_1$')
ylabel('$x_2$')
box on

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
plot(tmp.mptPolytope.P, 'color', 'cyan');
hold on
end

box on
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$\mathcal{U}([x])$')

%%
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 5 3]);
saveas(gcf,'img/e.eps','epsc')
saveas(gcf,'img/e.png')

%% Area
ar_true = 0.05*2*0.01*2;
ar = 0;
for i = 1:size(P,1)
ar = ar + Polyhedron(P{i}).volume;
end
ar;
ar2 = 0;
for i = 1:size(S,1)
ar2 = ar2 + S{i}.mptPolytope.P.volume;
end

ar/ar_true
ar2/ar_true
