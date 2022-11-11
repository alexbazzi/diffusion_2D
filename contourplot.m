function contourplot(T, colorMode, sc1, sc2)
figure('Name','Temperature Contour','NumberTitle','off');
[h, h] = contourf(flipud(T), 999);
if ~strcmp(colorMode, 'none')
    colormap(flipud(colorMode));
end
if ~strcmp(sc1, 'none') || ~strcmp(sc2, 'none')
    caxis([sc1;sc2]);
end
set(h,'LineStyle','none');
end