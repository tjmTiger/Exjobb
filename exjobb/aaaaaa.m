% format Watts Strogatz avergae degree figure from DF
f = gcf;

% Trivial Sol.
% Remove extra entries on Watts Strogatz
a = f.Children(2).Children(1).YData;
b = f.Children(2).Children(1).XData;
f.Children(2).Children(1).YData = [a(1) a(3) a(5)];
f.Children(2).Children(1).XData = [b(1) b(3) b(5)];

% remove triv sol for output feedback
f.Children(2).Children(2).YData = [];
f.Children(2).Children(2).YData = [];

% Runtime & cost
% Remove extra entries on Watts Strogatz
for i = 3:4
    for j = 1:2
        f.Children(i).Children(j).YData
        a = f.Children(i).Children(j).YData;
        b = f.Children(i).Children(j).XData;
        f.Children(i).Children(j).YData = [a(1) a(3) a(5)];
        f.Children(i).Children(j).XData = [b(1) b(3) b(5)];
    end
end