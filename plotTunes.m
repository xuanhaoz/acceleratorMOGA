function plotTunes(Q)
    % simple wrapper to plot tune footprint array
    figure(11021)
    plot(Q(:,1),Q(:,2),'Marker','+','LineWidth',2);
    grid on
end