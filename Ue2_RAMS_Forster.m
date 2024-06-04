clear;

function ecg_analysis(file_path, intervals)
    % Load the ECG data
    data = load(file_path);
    ecg_signal = data.ecg1;

    % Initialize results
    results = [];

    % Analyze each interval
    for i = 1:size(intervals, 1)
        start_time = intervals(i, 1) + 1; % Convert to 1-based indexing
        end_time = intervals(i, 2);

        if start_time > 0 && end_time <= length(ecg_signal)
            segment = ecg_signal(start_time:end_time);

            % Detect R peaks using a simple threshold method
            r_peaks = find_r_peaks(segment, 0.5);

            if length(r_peaks) > 1
                % Calculate HRV measures
                nn_intervals = diff(r_peaks);
                sdnn = std(nn_intervals);
                sdsd = std(diff(nn_intervals));

                % Calculate HR
                hr = 60000 ./ nn_intervals;
                avg_hr = mean(hr);
                min_hr = min(hr);
                max_hr = max(hr);

                % Store results
                results = [results; avg_hr, min_hr, max_hr, sdnn, sdsd];
            else
                % Store warning if not enough data
                results = [results; NaN, NaN, NaN, NaN, NaN];
            end
        else
            % Store warning if interval is out of bounds
            results = [results; NaN, NaN, NaN, NaN, NaN];
        end
    end

    % Summarize results
    summarize_results(results);
end

function r_peaks = find_r_peaks(segment, threshold)
    refractory_period = 200; % Minimum time between peaks (in samples)
    peaks = find(segment > threshold);
    if isempty(peaks)
        r_peaks = [];
    else
        r_peaks = peaks(1);
        for i = 2:length(peaks)
            if peaks(i) - r_peaks(end) > refractory_period
                r_peaks = [r_peaks; peaks(i)];
            end
        end
    end
end

function summarize_results(results)
    % Check if more than 7 intervals were analyzed
    valid_results = results(~isnan(results(:, 1)), :);
    if size(valid_results, 1) > 7
        avg_hr_all = mean(valid_results(:, 1));
        min_hr_all = min(valid_results(:, 2));
        max_hr_all = max(valid_results(:, 3));
        mean_sdnn = mean(valid_results(:, 4));
        std_sdnn = std(valid_results(:, 4));
        mean_sdsd = mean(valid_results(:, 5));
        std_sdsd = std(valid_results(:, 5));

        fprintf('Overall Statistics:\n');
        fprintf('Average HR: %.2f bpm\n', avg_hr_all);
        fprintf('HR Min: %.2f bpm, HR Max: %.2f bpm\n', min_hr_all, max_hr_all);
        fprintf('SDNN: Mean = %.2f ms, Std = %.2f ms\n', mean_sdnn, std_sdnn);
        fprintf('SDSD: Mean = %.2f ms, Std = %.2f ms\n', mean_sdsd, std_sdsd);
    else
        for i = 1:size(results, 1)
            if isnan(results(i, 1))
                fprintf('Interval %d: Not enough data to calculate HRV measures.\n', i);
            else
                fprintf('Interval %d: Avg HR = %.2f bpm, Min HR = %.2f bpm, Max HR = %.2f bpm, SDNN = %.2f ms, SDSD = %.2f ms\n', ...
                    i, results(i, 1), results(i, 2), results(i, 3), results(i, 4), results(i, 5));
            end
        end
    end
end

% Example usage with the provided data
file_path = 'ecg1.mat';
intervals = [0 10000; 10000 20000; 20000 30000; 30000 40000; 40000 50000; 50000 60000; 60000 70000; 70000 80000];
ecg_analysis(file_path, intervals);


