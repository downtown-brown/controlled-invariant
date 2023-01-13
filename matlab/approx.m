A = [1 1; 0 0.9]

for i = 10:100
    system = LTISystem('A', A - 10^(-i/20), 'B', [1; 0.5]);
    system.x.min = [-5; -5];
    system.x.max = [5; 5];
    system.u.min = -1;
    system.u.max = 1;
    InvSet = system.invariantSet();
    plot(InvSet, 'color', [0 .5 1/log(i-7)])
    hold on
end
