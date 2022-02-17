%% Plot Convex Hulls With Time
clear all;
close all;
clc;

for kk = 1:10
    figure();
    hold on;
    
    majorR = false;
    
    for hh = 1:2
        if hh == 1
            C_prev = [randi(10)*rand(), randi(10)*rand()];
            CP = C_prev;
        else
            startIndex = 1;
            while startIndex == 1 || startIndex == 6
                startIndex = randi(6);
            end
        
            C_prev = CP(startIndex, :);
            CP = C_prev;
        end
        
        R = randi(10)*rand();
        for ii = 1:6
            R_new = randi(10)*rand();
            
            if majorR
                while R_new > R
                    R_new = randi(10)*rand();
                end
            end
            
            R = R_new;
            
            C_angle = randi(360)*2*pi/360;
            C_coord = C_prev + [R*cos(C_angle), R*sin(C_angle)];
            plot(C_prev(1), C_prev(2), 'o', 'LineWidth', 1.5);
        
            if ii < 6
                plot([C_prev(1), C_coord(1)], [C_prev(2), C_coord(2)], '--');
                CP(ii + 1, :) = C_coord';
            else
                plot([C_prev(1), CP(1,1)], [C_prev(2), CP(1,2)], '--');
            end
        
            C_prev = C_coord;
        end
        
        if hh == 1
            CP_drone = CP;
        else
            CP_obs = CP;
        end
    end

    % Plot Curves %%%
    bin = [];
    r = [];

    for k = 0:1:5
        bin(k+1) = factorial(5)/(factorial(5-k)*factorial(k));
    end
    
    for tt = 0:0.001:1
        for k = 1:1:6
            b(k) = bin(k)*(1-tt)^(5-(k-1))*tt^(k-1);
        end
        r = cat(1, r, b);
    end
    Rb_drone = r*CP_drone;
    Rb_obs = r*CP_obs;

    plot(Rb_drone(:,1), Rb_drone(:,2), 'LineWidth', 1.5);
    plot(Rb_obs(:,1), Rb_obs(:,2), 'LineWidth', 1.5);
    %%%%

    % Sample Time Vector %%%
    for hh = 1:2
        CP_u(1) = 0;
        for ii = 2:5
            actual_u = rand();
            while actual_u < CP_u(ii-1)
                actual_u = rand();
            end
            CP_u(ii) = actual_u;
        end
        CP_u(6) = 1;

        if hh == 1
            CP_drone_u = CP_u';
        else
            CP_obs_u = CP_u';
        end
    end
    %%%%

    % Check Collision %%%
    bin = [];
    r = [];

    for k = 0:1:5
        bin(k+1) = factorial(5)/(factorial(5-k)*factorial(k));
    end
    
    for tt = 0:0.001:1
        for k = 1:1:6
            bu(k) = bin(k)*(1-tt)^(5-(k-1))*tt^(k-1);
        end

        uu_drone = bu*CP_drone_u;
        uu_obs = bu*CP_obs_u;

        for k = 1:1:6
            b_drone(k) = bin(k)*(1-uu_drone)^(5-(k-1))*uu_drone^(k-1);
            b_obs(k) = bin(k)*(1-uu_obs)^(5-(k-1))*uu_obs^(k-1);
        end

        p_drone = b_drone*CP_drone;
        p_obs = b_obs*CP_obs;

        if norm(p_drone - p_obs) < 0.5
            sprintf('Collision in curves %d', kk)
            break;
        end
    end
    hold off;
    %%%%

    CP_drone_mix = cat(2, CP_drone, CP_drone_u);
    CP_obs_mix = cat(2, CP_obs, CP_obs_u);

    [k_drone, av_drone] = convhull(CP_drone_mix);
    [k_obs, av_obs] = convhull(CP_obs_mix);

    figure();
    trisurf(k_drone, CP_drone_mix(:,1), CP_drone_mix(:,2), CP_drone_mix(:,3), 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    hold on;
    trisurf(k_obs, CP_obs_mix(:,1), CP_obs_mix(:,2), CP_obs_mix(:,3), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold off;
end
