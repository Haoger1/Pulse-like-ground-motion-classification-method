function results = classifyPulselikeGM(SWresults, displayFlag)
%CLASSIFYPULSELIKEGM Identify and classify pulse-like ground motions.
%
% Description:
%   This function identifies and classifies pulse-like ground motions based
%   on decomposed shock-waveform (SW) components. Each record is classified
%   as SPGM, DPGM, MPGM, or NPGM using peak-related parameters, duration-related
%   parameters, principal-component-based criteria, energy proportion (EP), PGV ratio,
%   and time-interval criteria.
%
% Syntax:
%   results = classifyPulseLikeGM(SWresults)
%   results = classifyPulseLikeGM(SWresults, displayFlag)
%
% Input:
%   SWresults
%       A cell array. The first column of the cell should contain the decomposed
%       results obtained from SWD. Each row of first column corresponds to each record. 
%
%   displayFlag
%       Optional logical value. If true, classification messages are displayed
%       in the MATLAB Command Window. Default is true.
%
% Required columns in each data matrix:
%   Column 1   : Time column of the original time series.
%   Column 2   : Velocity column of the original time series.
%   Column 3   : Energy proportion of each decomposed SW component relative
%                to the original signal. Each row corresponds to one decomposed
%                SW component. Therefore, the number of valid rows in Column 3
%                equals the number of decomposed SW components.
%   Column 3+j : The j-th decomposed SW component, where j = 1, 2, ..., n.
%                Therefore, Column 4 stores the first decomposed SW component,
%                Column 5 stores the second decomposed SW component, and so on.
%
% Output:
%   results
%       A structure array containing the classification result for each record.
%
%       results(i).recordIndex
%       results(i).class
%       results(i).subclass
%       results(i).selectedComponents
%       results(i).firstEnergyProportion
%       results(i).message
%
%   The classification categories are: 
%       SPGM: single pulse-like ground motion
%       DPGM: double pulse-like ground motion
%       MPGM: multiple pulse-like ground motion with more than two pulses (denoted here only since MPGM includes DPGM in original paper)
%       NPGM: non-pulse-like ground motion
% Notes:
%   - The function findpeaks requires MATLAB Signal Processing Toolbox.
%   - The algorithm removes rows containing zeros or NaN values in Column 3
%     to estimate the number of decomposed SW components.

    if nargin < 2
        displayFlag = true;
    end

    nRecords = length(SWresults);

    results = struct( ...
        'recordIndex', cell(nRecords, 1), ...
        'class', cell(nRecords, 1), ...
        'subclass', cell(nRecords, 1), ...
        'selectedComponents', cell(nRecords, 1), ...
        'firstEnergyProportion', cell(nRecords, 1), ...
        'message', cell(nRecords, 1));

    for i = 1:nRecords

        data = SWresults{i, 1};

        results(i).recordIndex = i;
        results(i).class = 'Unclassified';
        results(i).subclass = '';
        results(i).selectedComponents = [];
        results(i).firstEnergyProportion = data(1, 3);
        results(i).message = '';

        % Basic data check.
        if size(data, 2) < 4
            msg = ['Record ', num2str(i), ': insufficient number of columns; classified as NPGM.'];
            results(i).class = 'NPGM';
            results(i).message = msg;
            if displayFlag
                disp(msg);
            end
            continue;
        end

        % Estimate the number of decomposed SW components.
        valid_energy_rows = find(data(:, 3) ~= 0 & ~isnan(data(:, 3)));
        n = length(valid_energy_rows);

        if n < 1
            msg = ['Record ', num2str(i), ': no valid decomposed SW component is found; classified as NPGM.'];
            results(i).class = 'NPGM';
            results(i).message = msg;
            if displayFlag
                disp(msg);
            end
            continue;
        end

        % Check PGV to avoid division by zero.
        PGV_original = max(abs(data(:, 2)));

        if PGV_original == 0 || isnan(PGV_original)
            msg = ['Record ', num2str(i), ': invalid original velocity signal; classified as NPGM.'];
            results(i).class = 'NPGM';
            results(i).message = msg;
            if displayFlag
                disp(msg);
            end
            continue;
        end

        % Normalize time and velocity.
        t_normal = data(:, 1) / max(data(:, 1));
        y_normal = data(:, 2) / PGV_original;

        % Normalize the first decomposed SW component.
        decomsignal_1_normal = data(:, 4) / PGV_original;

        % -----------------------------------------------------------------
        % Case 1: EP of the first decomposed SW component is greater than 0.7
        % -----------------------------------------------------------------
        if data(1, 3) > 0.7

            [isValid, NP_log, PD] = getComponentFeatures(decomsignal_1_normal, t_normal);

            if ~isValid
                msg = ['Record ', num2str(i), ...
                       ': invalid first decomposed SW component; classified as NPGM.'];

                results(i).class = 'NPGM';
                results(i).message = msg;

            else
                x = PD + 0.642 * NP_log;

                if x < 1.106
                    msg = ['Record ', num2str(i), ...
                           ': the energy proportion of the first decomposed SW component is > 0.7; classified as SPGM.'];

                    results(i).class = 'SPGM';
                    results(i).subclass = 'High-energy first component';
                    results(i).selectedComponents = 1;
                    results(i).message = msg;
                else
                    msg = ['Record ', num2str(i), ...
                           ': the energy proportion of the first decomposed SW component is > 0.7; classified as NPGM.'];

                    results(i).class = 'NPGM';
                    results(i).message = msg;
                end
            end

            if displayFlag
                disp(msg);
            end

        % -----------------------------------------------------------------
        % Case 2: EP of the first decomposed SW component is between 0.3 and 0.7
        % -----------------------------------------------------------------
        elseif data(1, 3) >= 0.3 && data(1, 3) <= 0.7

            [isValid1, NP_log_1, PD_1] = getComponentFeatures(decomsignal_1_normal, t_normal);

            if ~isValid1
                msg = ['Record ', num2str(i), ...
                       ': invalid first decomposed SW component; classified as NPGM.'];

                results(i).class = 'NPGM';
                results(i).message = msg;

                if displayFlag
                    disp(msg);
                end

                continue;
            end

            PC1_1 = 0.968 * NP_log_1 + 0.251 * PD_1;
            PC1_pie_1 = 0.94 * NP_log_1 + 0.341 * PD_1;

            PGV_resi_normal = max(abs(y_normal - decomsignal_1_normal));
            PGVR = PGV_resi_normal / max(abs(y_normal));

            x = PGVR + 1.925e-11 * exp(PC1_1 / 0.0344);

            if x < 0.748
                msg = ['Record ', num2str(i), ...
                       ': the energy proportion of the first decomposed SW component is between 0.3 and 0.7; classified as SPGM.'];

                results(i).class = 'SPGM';
                results(i).subclass = 'SPGM_1';
                results(i).selectedComponents = 1;
                results(i).message = msg;

                if displayFlag
                    disp(msg);
                end

            else
                % If the second decomposed component does not exist, classify as NPGM.
                if n < 2 || size(data, 2) < 5
                    msg = ['Record ', num2str(i), ...
                           ': fewer than two decomposed SW components are available; classified as NPGM.'];

                    results(i).class = 'NPGM';
                    results(i).message = msg;

                    if displayFlag
                        disp(msg);
                    end

                    continue;
                end

                decomsignal_2_normal = data(:, 5) / PGV_original;
                [isValid2, NP_log_2, PD_2] = getComponentFeatures(decomsignal_2_normal, t_normal);

                if ~isValid2
                    msg = ['Record ', num2str(i), ...
                           ': invalid second decomposed SW component; classified as NPGM.'];

                    results(i).class = 'NPGM';
                    results(i).message = msg;

                    if displayFlag
                        disp(msg);
                    end

                    continue;
                end

                PC1_pie_2 = 0.94 * NP_log_2 + 0.341 * PD_2;

                x = PC1_pie_2 - 7.279e-5 * exp(PC1_pie_1 / 0.0753);

                if x > 1.004
                    msg = ['Record ', num2str(i), ...
                           ': the energy proportion of the first decomposed SW component is between 0.3 and 0.7; classified as SPGM.'];

                    results(i).class = 'SPGM';
                    results(i).subclass = 'SPGM_2';
                    results(i).selectedComponents = 2;
                    results(i).message = msg;

                    if displayFlag
                        disp(msg);
                    end

                else
                    xx = PC1_pie_2 + 0.0115 * exp(PC1_pie_1 / 0.313);

                    if xx < 1.383
                        [classLabel, subclassLabel, selectedComponents, msg] = ...
                            classifySelectedComponents(data, n, PGV_original, t_normal, y_normal, i, ...
                            'between 0.3 and 0.7');

                        results(i).class = classLabel;
                        results(i).subclass = subclassLabel;
                        results(i).selectedComponents = selectedComponents;
                        results(i).message = msg;

                    else
                        msg = ['Record ', num2str(i), ...
                               ': the energy proportion of the first decomposed SW component is between 0.3 and 0.7; classified as NPGM.'];

                        results(i).class = 'NPGM';
                        results(i).message = msg;
                    end

                    if displayFlag
                        disp(msg);
                    end
                end
            end

        % -----------------------------------------------------------------
        % Case 3: EP of the first decomposed SW component is less than 0.3
        % -----------------------------------------------------------------
        else

            [classLabel, subclassLabel, selectedComponents, msg] = ...
                classifySelectedComponents(data, n, PGV_original, t_normal, y_normal, i, '< 0.3');

            results(i).class = classLabel;
            results(i).subclass = subclassLabel;
            results(i).selectedComponents = selectedComponents;
            results(i).message = msg;

            if displayFlag
                disp(msg);
            end
        end
    end
end





function [isValid, NP_log, PD] = getComponentFeatures(component_normal, t_normal)
%GETCOMPONENTFEATURES Calculate NP_log and normalized pulse duration.

    isValid = false;
    NP_log = NaN;
    PD = NaN;

    [maxima, ~] = findpeaks(abs(component_normal), 'MinPeakHeight', 0.02);
    num_maxima = length(maxima);

    if num_maxima == 0
        return;
    end

    active_indices = find(abs(component_normal) >= 0.01);

    if isempty(active_indices)
        return;
    end

    NP_log = log10(num_maxima);

    first_active_t = t_normal(active_indices(1));
    last_active_t = t_normal(active_indices(end));

    PD = last_active_t - first_active_t;

    isValid = true;
end





function [classLabel, subclassLabel, selectedComponents, msg] = ...
    classifySelectedComponents(data, n, PGV_original, t_normal, y_normal, recordIndex, caseText)
%CLASSIFYSELECTEDCOMPONENTS Classify records based on selected SW components.

    pass_count = 0;
    pass_indices = [];

    for j = 1:n

        component_col = 3 + j;

        if component_col > size(data, 2)
            continue;
        end

        decomposed_signal_normal = data(:, component_col) / PGV_original;

        [isValid, NP_log, PD] = getComponentFeatures(decomposed_signal_normal, t_normal);

        if ~isValid
            continue;
        end

        PGV_PGVori = max(abs(decomposed_signal_normal)) / max(abs(y_normal));

        x = PD + 1.21 * (NP_log - 0.845)^2;

        if x < 0.6 && PGV_PGVori > 0.5
            pass_count = pass_count + 1;
            pass_indices = [pass_indices; j];
        end
    end

    selectedComponents = pass_indices;

    if pass_count == 1

        comp1 = data(:, 3 + pass_indices(1)) / PGV_original;

        PGV_resi_normal = max(abs(y_normal - comp1));
        PGVR = PGV_resi_normal / max(abs(y_normal));

        [isValid, NP_log, PD] = getComponentFeatures(comp1, t_normal);

        if ~isValid
            classLabel = 'NPGM';
            subclassLabel = '';
            selectedComponents = [];
            msg = ['Record ', num2str(recordIndex), ...
                   ': the energy proportion of the first decomposed SW component is ', caseText, ...
                   '; invalid selected component; classified as NPGM.'];
            return;
        end

        PC1_piepie = 0.941 * NP_log + 0.338 * PD;

        x = PGVR + 3.16e-3 * exp(PC1_piepie / 0.178);

        if x < 0.678
            classLabel = 'SPGM';
            subclassLabel = 'Selected single component';

            msg = ['Record ', num2str(recordIndex), ...
                   ': the energy proportion of the first decomposed SW component is ', caseText, ...
                   '; n = 1 and component ', num2str(pass_indices(1)), ...
                   ' is selected; classified as SPGM.'];
        else
            classLabel = 'NPGM';
            subclassLabel = '';
            selectedComponents = [];

            msg = ['Record ', num2str(recordIndex), ...
                   ': the energy proportion of the first decomposed SW component is ', caseText, ...
                   '; n = 1 but classified as NPGM.'];
        end

    elseif pass_count == 2

        comp1 = data(:, 3 + pass_indices(1)) / PGV_original;
        comp2 = data(:, 3 + pass_indices(2)) / PGV_original;

        PGV_resi_normal = max(abs(y_normal - comp1 - comp2));
        PGVR = PGV_resi_normal / max(abs(y_normal));

        [~, idx_y1] = max(abs(comp1));
        [~, idx_y2] = max(abs(comp2));

        TI = abs(t_normal(idx_y1) - t_normal(idx_y2));

        x = PGVR + 1.47 * exp(-TI / 0.0229);

        if x < 0.619
            classLabel = 'MPGM';
            subclassLabel = 'DPGM';

            msg = ['Record ', num2str(recordIndex), ...
                   ': the energy proportion of the first decomposed SW component is ', caseText, ...
                   '; n = 2 and components ', num2str(pass_indices(1)), ...
                   ' and ', num2str(pass_indices(2)), ...
                   ' are selected; classified as MPGM (DPGM).'];
        else
            classLabel = 'NPGM';
            subclassLabel = '';
            selectedComponents = [];

            msg = ['Record ', num2str(recordIndex), ...
                   ': the energy proportion of the first decomposed SW component is ', caseText, ...
                   '; n = 2 but classified as NPGM.'];
        end

    elseif pass_count > 2

        decompsignal_sum_normal = zeros(size(data, 1), 1);

        for k = 1:length(pass_indices)
            decompsignal_sum_normal = decompsignal_sum_normal ...
                + data(:, 3 + pass_indices(k)) / PGV_original;
        end

        PGV_resi_normal = max(abs(y_normal - decompsignal_sum_normal));
        PGVR = PGV_resi_normal / max(abs(y_normal));

        TI_all = [];

        for k = 1:length(pass_indices) - 1

            comp_current = data(:, 3 + pass_indices(k)) / PGV_original;
            comp_next = data(:, 3 + pass_indices(k + 1)) / PGV_original;

            [~, idx_y1] = max(abs(comp_current));
            [~, idx_y2] = max(abs(comp_next));

            TI = abs(t_normal(idx_y1) - t_normal(idx_y2));
            TI_all = [TI_all; TI];
        end

        pass_indices_str = sprintf('%d, ', pass_indices);
        pass_indices_str = pass_indices_str(1:end-2);

        if all(TI_all > 0.05) && PGVR < 0.6
            classLabel = 'MPGM';
            subclassLabel = 'Multiple selected components';

            msg = ['Record ', num2str(recordIndex), ...
                   ': the energy proportion of the first decomposed SW component is ', caseText, ...
                   '; n = ', num2str(length(pass_indices)), ...
                   ' and pass_indices = ', pass_indices_str, ...
                   '; classified as MPGM.'];
        else
            classLabel = 'NPGM';
            subclassLabel = '';
            selectedComponents = [];

            msg = ['Record ', num2str(recordIndex), ...
                   ': the energy proportion of the first decomposed SW component is ', caseText, ...
                   '; n = ', num2str(length(pass_indices)), ...
                   ' but classified as NPGM.'];
        end

    else

        classLabel = 'NPGM';
        subclassLabel = '';
        selectedComponents = [];

        msg = ['Record ', num2str(recordIndex), ...
               ': the energy proportion of the first decomposed SW component is ', caseText, ...
               '; n = 0, therefore classified as NPGM.'];
    end
end

