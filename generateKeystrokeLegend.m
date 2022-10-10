function generateKeystrokeLegend(brushing)

fig = findall(0, 'Type', 'figure', 'name', 'Keystroke Command Legend');
if isempty(fig)
    fig = figure('Name', 'Keystroke Command Legend');
end
fig.Position = brushing.params.commandLegendPos;


keystrokes = {'0', '1-8', 'd', 'z', 'x', 'a', 'u'};
functs = {'remove label', 'label as #1-8', 'delete', 'zoom ON', 'zoom OFF', 'associate', 'undo'};

vertLines = [0, 2, 8];
horzLines = 1:numel(keystrokes);

xdata = [1, 5];
ydata = (1:numel(keystrokes)) + .5;

axis ij
hold on
for i = 1:length(vertLines)
    plot([vertLines(i),vertLines(i)], [1, 8], 'k')
end
for i = 1:length(horzLines)
    plot([0, 8], [horzLines(i),horzLines(i)], 'k')
end
for i = 1:length(keystrokes)
    text(xdata(1), ydata(i), keystrokes{i}, 'HorizontalAlignment', 'center', 'FontSize', 13)
end
for i = 1:length(functs)
    text(xdata(2), ydata(i), functs{i}, 'HorizontalAlignment', 'center', 'FontSize', 11)
end
fig.Children.XTick = [];
fig.Children.YTick = [];

xlim([0, 8]);
ylim([1, 8])
hold off


