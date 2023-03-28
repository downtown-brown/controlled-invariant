figure
hold on
% loops = 14;
% M(loops) = struct('cdata',[],'colormap',[]);
for i = 1:size(S, 2)
    for s = 1:size(S{i}, 1)
        plot(S{i}{s},[1 2],'FaceColor',[.1 .1 .6],'Filled',true,'EdgeColor','black');
    end
    
    for s = 1:size(N{i}, 1)
        plot(N{i}{s},[1 2],'FaceColor',[.6 .1 .1],'Filled',true,'EdgeColor','black');
    end
    
    for s = 1:size(E{i}, 1)
        plot(E{i}{s},[1 2],'FaceColor',[.1 .6 .1],'Filled',true,'EdgeColor','black');
    end
%     if i > 1
%     for s = 1:size(E{i-1}, 1)
%         plot(E{i-1}{s},[1 2],'FaceColor',[.6 .1 .1],'Filled',true,'EdgeColor','black');
%     end
%     end

%     M(i) = getframe;
keyboard
end