function drawObstacles()
    % Obstacles
    obs1 = [2,0,2,2];
    obs2 = [2,3,2,2];
    obs3 = [2,6,2,2];
    obs4 = [5,1,1,6];
    
    rectangle('Position', obs1, 'FaceColor', [1 1 1], 'EdgeColor','k', 'LineWidth', 3, 'Curvature',0.2)
    rectangle('Position', obs2, 'FaceColor', [1 1 1], 'EdgeColor','k', 'LineWidth', 3, 'Curvature',0.2)
    rectangle('Position', obs3, 'FaceColor', [1 1 1], 'EdgeColor','k', 'LineWidth', 3, 'Curvature',0.2)
    rectangle('Position', obs4, 'FaceColor', [1 1 1], 'EdgeColor','k', 'LineWidth', 3, 'Curvature',0.2)
end