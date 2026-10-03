function [maxh_matrix, AngleALL] = computeMaxh_nonlinear(clstack, corrstack, K, p, method, sigma, a)
% Nonlinear weighted voting for computing common line reliability
% Implements: w_ij^k = f(H_ik) * f(H_jk)  (product of two supporting CL reliabilities)
% method: 'sigmoid' | 'tanh' | 'linear' | 'power' | 'gaussian'
% p    : percentile for threshold (applied to maxh_matrix_old)
% a    : steepness for sigmoid/tanh, exponent for power
% sigma: width for gaussian

% Initialize
maxh_matrix = zeros(K, K);
AngleALL    = zeros(K, K);
T   = 60;
tho = 180 / T;
count = 0;

% Stage 1: uniform-weight initialization
maxh_matrix_old = computeMaxhandSts_ang(clstack, corrstack, K);

% Iteration control
max_iter               = 30;
patience               = 5;
relative_std_threshold = 0.01;
error_history          = [];

% Threshold fixed from first iteration (percentile of maxh)
Thresh = prctile(maxh_matrix_old(:), p);
fprintf('Thresh = %.4f  (p=%d percentile of maxh_matrix_old)\n', Thresh, p);

sum_countW_old     = 0;
sum_maxh_matrix_old = sum(maxh_matrix_old(:));

while count < max_iter
    sum_countW = 0;

    for k1 = 1:K-1
        for k2 = k1+1:K
            h      = zeros(1, T);
            countW = 0;
            sum_w  = 0;

            for k3 = 1:K
                if k3 == k1 || k3 == k2
                    continue;
                end

                % Gate: both supporting common lines must be reliable
                if maxh_matrix_old(k1,k3) > Thresh && maxh_matrix_old(k2,k3) > Thresh

                    % === Weight: f(H_ik) * f(H_jk)  (paper Eq.7-8) ===
                    H_ik = maxh_matrix_old(k1, k3);
                    H_jk = maxh_matrix_old(k2, k3);

                    if strcmp(method, 'sigmoid')
                        w_ik = 1 / (1 + exp(-a * (H_ik - Thresh)));
                        w_jk = 1 / (1 + exp(-a * (H_jk - Thresh)));

                    elseif strcmp(method, 'tanh')
                        w_ik = 0.5 * (1 + tanh(a * (H_ik - Thresh)));
                        w_jk = 0.5 * (1 + tanh(a * (H_jk - Thresh)));

                    elseif strcmp(method, 'linear')
                        % linear: w = H directly (maxh is already in [0,1])
                        w_ik = H_ik;
                        w_jk = H_jk;

                    elseif strcmp(method, 'power')
                        norm_ik = (H_ik - Thresh) / (1 - Thresh);
                        norm_jk = (H_jk - Thresh) / (1 - Thresh);
                        w_ik = norm_ik ^ a;
                        w_jk = norm_jk ^ a;

                    elseif strcmp(method, 'gaussian')
                        w_ik = exp(-(H_ik - Thresh)^2 / (2 * sigma^2));
                        w_jk = exp(-(H_jk - Thresh)^2 / (2 * sigma^2));

                    else
                        error('Unknown method: %s', method);
                    end

                    wk = w_ik * w_jk;   % product (paper Eq.8)
                    countW = countW + 1;
                    sum_w  = sum_w  + wk;

                else
                    wk = 0;
                end

                % Compute dihedral angle estimate
                [angle12, flag] = ComputeAngleCorr( ...
                    clstack(k1,k2), clstack(k2,k1), corrstack(k1,k2), ...
                    clstack(k1,k3), clstack(k3,k1), corrstack(k1,k3), ...
                    clstack(k2,k3), clstack(k3,k2), corrstack(k2,k3));

                if flag && wk > 0
                    ang = (1:T) * 180.0 / T;
                    h   = h + wk * exp(-(ang - angle12).^2 / (2 * tho * tho));
                end
            end

            % Normalize by total weight (paper Eq.5)
            if countW == 0 || sum_w == 0
                maxh_matrix(k1,k2) = 0;
                maxh_matrix(k2,k1) = 0;
                continue;
            else
                h = h ./ sum_w;
                sum_countW = sum_countW + countW;
            end

            [maxh, idx] = max(h);
            ang = idx * tho - tho / 2;

            maxh_matrix(k1,k2) = maxh;
            maxh_matrix(k2,k1) = maxh;
            AngleALL(k1,k2)    = ang;
            AngleALL(k2,k1)    = ang;
        end
    end

    % Convergence check
    current_error = max(max(abs(maxh_matrix - maxh_matrix_old)));
    error_history = [error_history; current_error];

    oscillation_detected = false;
    if length(error_history) >= patience
        recent_errors = error_history(end-patience+1:end);
        recent_mean   = mean(recent_errors);
        recent_std    = std(recent_errors);
        if recent_mean > 0 && (recent_std / recent_mean) < relative_std_threshold
            fprintf('Oscillation detected: mean=%.6f  rel_std=%.2f%%\n', ...
                recent_mean, (recent_std/recent_mean)*100);
            oscillation_detected = true;
        end
    end

    fprintf('[iter %d]  max_old=%.4f  max_new=%.4f  error=%.6f\n', ...
        count+1, max(maxh_matrix_old(:)), max(maxh_matrix(:)), current_error);

    % Convergence criterion (paper Eq.9): error < 1% of current max
    converged = current_error < 0.01 * max(maxh_matrix_old(:));

    if converged || oscillation_detected
        fprintf('Converged at iter %d\n', count+1);
        break;
    else
        sum_maxh_matrix = sum(maxh_matrix(:));
        fprintf('sum_countW: %d -> %d\n', sum_countW_old, sum_countW);
        fprintf('sum_maxh:   %.4f -> %.4f\n', sum_maxh_matrix_old, sum_maxh_matrix);
        maxh_matrix_old     = maxh_matrix;
        sum_maxh_matrix_old = sum_maxh_matrix;
        sum_countW_old      = sum_countW;
        count = count + 1;
    end
end

fprintf('========== Done: %d iterations ==========\n', count);
fprintf('Final error: %.6f\n', current_error);

end