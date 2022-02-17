%% Plot Convex Hulls With Time
clear all;
close all;
clc;

% Curve Orders
m = 3;
n = 3;

P = sym('p', [m+1, 1]);
U = sym('u', [n+1, 1]);

A = P;

for kk = 1:m
    for ii = 0:m-kk
        for jj = 0:kk*n
            iIdx = ii+1;
            kIdx = kk+1;
            jIdx = jj+1;

            sum = 0;
            for ll = max([0, jj-n]):min([jj, kk*n-n])
                lIdx = ll+1;

                w1 = binomial(kk*n-n, ll);
                w2 = binomial(n, jj-ll);
                w3 = ((1-U(jIdx-lIdx+1,1))*A(iIdx, kIdx-1, lIdx)) + ...
                     (U(jIdx-lIdx+1,1)*A(iIdx+1, kIdx-1, lIdx));

                sum = sum + (w1*w2*w3);
            end

            A(iIdx, kIdx, jIdx) = (1/binomial(kk*n, jj))*sum;
        end
    end
end

C(:,1) = A(1, m+1, 1:((m*n)+1));
C = simplify(C); % New Control Points

for kk = 1:2
    figure();
    hold on;
    
    majorR = false;
    
    for hh = 1:2
        if hh == 1
            C_prev = [randi(10)*rand(), randi(10)*rand()];
            CP = C_prev;
        else
            startIndex = 1;
            while startIndex == 1 || startIndex == m+1
                startIndex = randi(m+1);
            end
        
            C_prev = CP(startIndex, :);
            CP = C_prev;
        end
        
        R = randi(10)*rand();
        for ii = 1:m+1
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
        
            if ii < m+1
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
    r = [];
    for k = 0:1:m
        bin(k+1) = factorial(m)/(factorial(m-k)*factorial(k));
    end
    
    for tt = 0:0.001:1
        for k = 1:1:m+1
            b(k) = bin(k)*(1-tt)^(m-(k-1))*tt^(k-1);
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
        for ii = 2:n
            actual_u = rand();
            while actual_u < CP_u(ii-1)
                actual_u = rand();
            end
            CP_u(ii) = actual_u;
        end
        CP_u(n+1) = 1;

        if hh == 1
            CP_drone_u = CP_u';
        else
            CP_obs_u = CP_u';
        end
    end
    %%%%

    % Check Collision %%%
    for k = 0:1:m
        bin(k+1) = factorial(m)/(factorial(m-k)*factorial(k));
    end

    for k = 0:1:n
        bin3(k+1) = factorial(n)/(factorial(n-k)*factorial(k));
    end

    for k = 0:1:(m*n)
        bin15(k+1) = factorial((m*n))/(factorial((m*n)-k)*factorial(k));
    end

    for k = 0:1:(m*n*2)
        bin30(k+1) = factorial((m*n*2))/(factorial((m*n*2)-k)*factorial(k));
    end
    
    collision = false;
    idx = 1;
    idx_ = 1;
    for tt = 0:0.001:1
        for k = 1:1:n+1
            bu(k) = bin3(k)*(1-tt)^(n-(k-1))*tt^(k-1);
        end

        uu_drone = bu*CP_drone_u;
        uu_obs = bu*CP_obs_u;

        for k = 1:1:m+1
            b_drone(k) = bin(k)*(1-uu_drone)^(m-(k-1))*uu_drone^(k-1);
            b_obs(k) = bin(k)*(1-uu_obs)^(m-(k-1))*uu_obs^(k-1);
            bb(k) = bin(k)*(1-tt)^(m-(k-1))*tt^(k-1);
        end

        for k = 1:1:(m*n)+1
            bb(k) = bin15(k)*(1-tt)^((m*n)-(k-1))*tt^(k-1);
        end

        for k = 1:1:(m*n*2)+1
            bb_sq(k) = bin30(k)*(1-tt)^((m*n*2)-(k-1))*tt^(k-1);
        end

        p_drone = b_drone*CP_drone;
        p_obs = b_obs*CP_obs;

        if norm(p_drone - p_obs) < 0.5 && not(collision)
            sprintf('Collision in curves %d', kk)
            collision = true;
        end

        trajWithTime_drone(idx, :) = [p_drone, uu_drone];
        trajWithTime_obs(idx, :) = [p_obs, uu_obs];

        bu_vect(idx, :) = bu;
        b_drone_vect(idx, :) = b_drone;
        bb_vect(idx, :) = bb; 
        bb_sq_vect(idx, :) = bb_sq; 

        idx = idx + 1;
    end
    hold off;

    %%%%%%
    % Compute Coposition P/U
    % Drone
    for ii = 1:2
        for jj = 1:size(P,1)
            assignin('caller', ['p', num2str(jj)], CP_drone(jj,ii));
        end

        for jj = 1:size(U,1)
            assignin('caller', ['u', num2str(jj)], CP_drone_u(jj,1));
        end
        
        CP_drone_new(:,ii) = subs(C); % Evaluate Composed Control Points
    end

    % Obstacle
    for ii = 1:2
        for jj = 1:size(P,1)
            assignin('caller', ['p', num2str(jj)], CP_obs(jj,ii));
        end

        for jj = 1:size(U,1)
            assignin('caller', ['u', num2str(jj)], CP_obs_u(jj,1));
        end
        
        CP_obs_new(:,ii) = subs(C); % Evaluate Composed Control Points
    end

    % U Degree Elevation n -> n*m
    r = (n*m)-n;
    for ii = 1:(n*m)+1
        CP_drone_u_new(ii, 1) = 0;
        CP_obs_u_new(ii, 1) = 0;

        k = ii-1;
        for jj = max([0, k-r]):min([n, k])
            factor = binomial(r, k-jj)*binomial(n, jj)/binomial(n*m, k);
            CP_drone_u_new(ii, 1) = CP_drone_u_new(ii, 1) + factor*CP_drone_u(jj+1);
        end

        for jj = max([0, k-r]):min([n, k])
            factor = binomial(r, k-jj)*binomial(n, jj)/binomial(n*m, k);
            CP_obs_u_new(ii, 1) = CP_obs_u_new(ii, 1) + factor*CP_obs_u(jj+1);
        end
    end

    CP_drone_mix = double([CP_drone_new, CP_drone_u_new]);
    CP_obs_mix = double([CP_obs_new, CP_obs_u_new]);

    [k_drone, av_drone] = convhull(CP_drone_mix);
    [k_obs, av_obs] = convhull(CP_obs_mix);

    figure()
    trisurf(k_drone, CP_drone_mix(:,1), CP_drone_mix(:,2), CP_drone_mix(:,3), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on;
    trisurf(k_obs, CP_obs_mix(:,1), CP_obs_mix(:,2), CP_obs_mix(:,3), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot3(trajWithTime_drone(:, 1), trajWithTime_drone(:, 2), trajWithTime_drone(:, 3), 'LineWidth', 1.5);
    plot3(trajWithTime_obs(:, 1), trajWithTime_obs(:, 2), trajWithTime_obs(:, 3), ':', 'LineWidth', 1.5);
    hold off;

    % Test
    Ru = bb_vect*CP_drone_mix(:, 3);
    Rb = bb_vect*CP_drone_mix(:, 1:2);

    Ru_obs = bb_vect*CP_obs_mix(:, 3);
    Rb_obs = bb_vect*CP_obs_mix(:, 1:2);

    hold on;
    plot3(Rb(:,1), Rb(:,2), Ru, 'LineWidth', 1.5);
    plot3(Rb_obs(:,1), Rb_obs(:,2), Ru_obs, 'LineWidth', 1.5);
    hold off;

    % Testing Norm^2 of Difference
    CP_diff = CP_drone_mix-CP_obs_mix;

    for kk = 1:(n*m*2)+1
        CP_diff_sq(kk,:) = [0, 0, 0];
        
        k = kk - 1;
        for jj = max([0, k-(n*m)]):min([n*m, k])
            factor = binomial(n*m, jj)*binomial(n*m, k-jj)/binomial(n*m*2, k);
            CP_diff_sq(kk,:) = CP_diff_sq(kk,:) + factor*(CP_diff(jj+1,:).*CP_diff(kk-jj,:));
        end
    end

    CP_norm_sq = CP_diff_sq(:, 1) + CP_diff_sq(:, 2) + CP_diff_sq(:, 3);
    R_sq = bb_sq_vect*CP_norm_sq;

    figure();
    hold on;
    plot([0:0.001:1], R_sq, 'LineWidth', 1.5);
    plot([0:0.001:1], (0.5^2)*ones(size([0:0.001:1])), 'LineWidth', 1.5);
%    ylim([0, 2]);
    hold off;

    figure();
    hold on;
    plot([1:size(CP_norm_sq)], CP_norm_sq, 'x', 'LineWidth', 1.5);
    plot([1:size(CP_norm_sq)], (0.5^2)*ones(size([1:size(CP_norm_sq)])), 'LineWidth', 1.5);
%     ylim([-2, 2]);
    hold off;
end

%% External Functions
function b = binomial(x, y)
    b = factorial(x)/(factorial(y)*factorial(x-y));
end
