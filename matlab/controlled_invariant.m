n = 2;
A = [0 1; .6 0];
r = [0.05; 0];

X = Interval([-10;-10],[10;10]);

u_h = [2,2]
u_l = [-2,-2]
U = Interval(u_l, u_h);

w_l = [-2;-1]/1
w_h = [2;1]/1
W = Interval(w_l,w_h);

B = [1;0];
B = eye(n);

T = 10;
x0 = [.1,-.3];

x = zeros(n,T);
x(:,1) = x0;
figure(1000)
clf
xlim([-11,11])
ylim([-11,11])
for t = 1:T-1
    u = rand(n,1)*2 - 1
    cvx_begin sdp
        variable x_l(n,1)
        variable x_h(n,1)
        maximize sum(x_h - x_l)
        subject to
            x_h >= x_l
            x_l >= [-10;-10]
            x_h <= [10;10]
            A*x_h + B*u + w_h <= x_h
            A*x_l + B*u + w_l >= x_l
    cvx_end
    disp([x_l, x_h])
    if strcmp(cvx_status,'Solved')
        plot_rect(x_l, x_h)
        hold on
    end
end
hold off

function plot_rect(x_l, x_h)
    w = x_h(1) - x_l(1);
    h = x_h(2) - x_l(2);
    pos = [x_l(1) x_l(2) w h];
    rectangle('Position', pos)
end