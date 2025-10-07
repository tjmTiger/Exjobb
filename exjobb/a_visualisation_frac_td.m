clear; clc;
mycolors = [
    238,64,53;
    243,155,54;
    123,192,67;
    3,146,207;
    17,0,255;
    175,56,255;
    ]./255;

n = 100;
figure;
c = 1;
for fD = 0.05:0.05:0.3
    C = [];
    for fT = 0.05:0.05:0.3
        T = fT*n;
        D = fD*n;
        V_in = T;
        C(end+1) = 2*V_in / (T+D);
    end
    hold on
    plot(0.05:0.05:0.3, C, "-o", 'DisplayName', num2str(fD), 'Color',mycolors(c,:))
    hold off
    c = c+1;
end

xlabel("FracTarg")
ylabel("Cost")
title("Worst case Cost$^{SF}$", "Interpreter",'latex')
lgd = legend('Location',"southeast");
title(lgd,'FracDist')

interpreter = 'latex';
all_text = findall(gcf, "-property", "Interpreter");
set(all_text, "Interpreter", interpreter)

all_text = findall(gcf, "-property", "TickLabelInterpreter");
set(all_text, "TickLabelInterpreter", interpreter)

all_text = findall(gcf, "-property", "Fontsize");
set(all_text, "Fontsize", 12)