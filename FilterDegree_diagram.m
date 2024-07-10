% Butterworth filter design
Amax=1;
Amax_ratio = linspace(1,4,1000);
ws_values = [1.1,1.2,1.3,1.4,1.5];
plot_n_vs_Amax_Amin_ratio(ws_values, Amin, Amax_ratio)

function plot_n_vs_Amax_Amin_ratio(ws_values, Amax, Amax_ratio)
    
    n = @(ws,Amin,Amax) ceil(log10( (10.^(Amin./10)-1)./(sqrt(10.^(Amax./10)-1).^2) )./(log10(ws)*2)); 

    % Calculate Amax based on the given Amin and Amax_ratio
    Amin = Amax_ratio * Amax;
    
    % Create a figure
    figure;
    hold on;

    % Loop through each ws value
    for ws = ws_values
        % Calculate n for the given Amin and Amax
        n_values = n(ws,Amin,Amax);

        % Plot n values against the ratio Amax/Amin
        plot(Amax_ratio, n_values, 'DisplayName', ['tran = ' num2str(ws-1)],'LineWidth',3);
    end

    % Add labels and title
    xlabel('Amin/Amax');
    ylabel('n');
    title('n vs Amax/Amin Ratio for different ws values');
    legend show;
    grid on;
    ylim([0, 11]);
    hold off;
end