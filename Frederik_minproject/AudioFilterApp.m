classdef AudioFilterApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure            matlab.ui.Figure
        LoadButton          matlab.ui.control.Button
        ApplyFiltersButton  matlab.ui.control.Button
        SaveButton          matlab.ui.control.Button
        LowGainSlider       matlab.ui.control.Slider
        MidGainSlider       matlab.ui.control.Slider
        HighGainSlider      matlab.ui.control.Slider
        WetnessSlider       matlab.ui.control.Slider
        OriginalAxes        matlab.ui.control.UIAxes
        LowGainLabel        matlab.ui.control.Label
        MidGainLabel        matlab.ui.control.Label
        HighGainLabel       matlab.ui.control.Label
        WetnessLabel        matlab.ui.control.Label
    end
    
    properties (Access = private)
        OriginalAudio % Original audio signal
        FilteredAudio % Filtered audio signal
        fs % Sampling frequency
    end
    
    methods (Access = private)
        
% Load audio file
function loadAudio(app)
    [file, path] = uigetfile('*.wav', 'Select an audio file');
    if isequal(file, 0)
        return;
    end
    [app.OriginalAudio, app.fs] = audioread(fullfile(path, file));
    
    % Clear the axes before plotting
    cla(app.OriginalAxes); % Clear the original audio axes
    applyFilters(app)
    plotPowerSpectrum(app, app.FilteredAudio, app.OriginalAxes, 'Filtered Audio');

end

% Apply filters
function applyFilters(app)
    disp('Applying filters...'); % Debug statement
    f_low = 100; % Cutoff frequency for low shelf filter
    f_mid = 1000; % Center frequency for peaking filter
    f_high = 8000; % Cutoff frequency for high shelf filter
    
    gain_low = app.LowGainSlider.Value;
    gain_mid = app.MidGainSlider.Value;
    gain_high = app.HighGainSlider.Value;
    
    wetness = app.WetnessSlider.Value / 100; % Convert to a scale from 0 to 1
    
    % Design filters
    disp('Designing filters...'); % Debug statement
    [b_low, a_low] = design_shelf_filter(app.fs, f_low, gain_low, 'low');
    [b_mid, a_mid] = design_peak_filter(app.fs, f_mid, gain_mid, 1);
    [b_high, a_high] = design_shelf_filter(app.fs, f_high, gain_high, 'high');
    
    % Apply filters in parallel
    disp('Filtering audio...'); % Debug statement
    y_low = filter(b_low, a_low, app.OriginalAudio);
    y_mid = filter(b_mid, a_mid, app.OriginalAudio);
    y_high = filter(b_high, a_high, app.OriginalAudio);
    
    % Combine filtered signals
    filtered = y_low + y_mid + y_high;
    
    % Normalize the combined signal to avoid clipping
    filtered = filtered / max(abs(filtered));
    
    % Apply wet/dry mix
    app.FilteredAudio = wetness * filtered + (1 - wetness) * app.OriginalAudio;
    
    % Normalize to avoid clipping
    app.FilteredAudio = app.FilteredAudio / max(abs(app.FilteredAudio));
    
    % Clear the axes before plotting the new power spectrum
    cla(app.OriginalAxes);
    
    

    % Plot the power spectrum of the filtered audio
    plotPowerSpectrum(app, app.FilteredAudio, app.OriginalAxes, 'Filtered Audio');

end










        
        % Save filtered audio
        function saveFilteredAudio(app)
        [file, path] = uiputfile('*.wav', 'Save Filtered Audio');
        if isequal(file, 0)
         return;
        end
    
    % Normalize the filtered audio to avoid clipping
    filtered_audio_normalized = app.FilteredAudio / max(abs(app.FilteredAudio));
    
    % Save the normalized filtered audio
    audiowrite(fullfile(path, file), filtered_audio_normalized, app.fs);
    uialert(app.UIFigure, 'Filtered audio saved successfully!', 'Success');
end
        
        % Plot power spectrum
function plotPowerSpectrum(app, signal, axesHandle, titleText)
    disp('Plotting power spectrum...'); % Debug statement
    
    % Calculate the power spectral density (PSD)
    [Pxx, F] = periodogram(signal, [], [], app.fs, 'power');
    
    % Convert power spectral density (PSD) to dB/Hz
    psd_dB_per_Hz = 10*log10(Pxx) - 10*log10(app.fs);
    
    % Smooth the power spectrum using a moving average filter
    windowSize = 10; % Choose an appropriate window size
    smooth_psd_dB_per_Hz = movmean(psd_dB_per_Hz, windowSize);
    
    % Plot the smoothed power spectrum in terms of frequency (Hz) by sound pressure level (dB/Hz)
    plot(axesHandle, F, smooth_psd_dB_per_Hz);
    
    title(axesHandle, titleText);
    xlabel(axesHandle, 'Frequency (Hz)');
    ylabel(axesHandle, 'Power Spectral Density (dB/Hz)');
    xlim(axesHandle, [0, app.fs/2]); % Set x-axis limit to Nyquist frequency
    grid(axesHandle, 'on');
end





    end

    % Callback methods
    methods (Access = private)

        % Button pushed function: LoadButton
        function LoadButtonPushed(app, event)
            loadAudio(app);
        end

        % Button pushed function: ApplyFiltersButton
        function ApplyFiltersButtonPushed(app, event)
            applyFilters(app);
        end

        % Button pushed function: SaveButton
        function SaveButtonPushed(app, event)
            saveFilteredAudio(app);
        end

        % Value changed function: LowGainSlider
        function LowGainSliderValueChanged(app, event)
        disp('LowGainSlider value changed.'); % Debug statement
        applyFilters(app);
        end

        % Value changed function: MidGainSlider
        function MidGainSliderValueChanged(app, event)
        disp('MidGainSlider value changed.'); % Debug statement
        applyFilters(app);
        end

        % Value changed function: HighGainSlider
        function HighGainSliderValueChanged(app, event)
        disp('HighGainSlider value changed.'); % Debug statement
        applyFilters(app);
        end

        % Value changed function: WetnessSlider
        function WetnessSliderValueChanged(app, event)
        disp('WetnessSlider value changed.'); % Debug statement
        applyFilters(app);
        end

    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create LoadButton
            app.LoadButton = uibutton(app.UIFigure, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.Position = [50 420 100 22];
            app.LoadButton.Text = 'Load Audio';

            % Create ApplyFiltersButton
            app.ApplyFiltersButton = uibutton(app.UIFigure, 'push');
            app.ApplyFiltersButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyFiltersButtonPushed, true);
            app.ApplyFiltersButton.Position = [50 380 100 22];
            app.ApplyFiltersButton.Text = 'Apply Filters';

            % Create SaveButton
            app.SaveButton = uibutton(app.UIFigure, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.Position = [50 340 100 22];
            app.SaveButton.Text = 'Save Audio';

            % Create LowGainSlider
            app.LowGainSlider = uislider(app.UIFigure);
            app.LowGainSlider.ValueChangedFcn = createCallbackFcn(app, @LowGainSliderValueChanged, true);
            app.LowGainSlider.Position = [200 190 200 3];
            app.LowGainSlider.Limits = [-20 20]; % Increase range for stronger effect
            app.LowGainSlider.Value = 0; % Set initial value
            app.LowGainSlider.MajorTicks = -20:2:20; % Set major ticks
            app.LowGainSlider.MinorTicks = -20:0.5:20; % Set minor ticks

            % Create LowGainLabel
            app.LowGainLabel = uilabel(app.UIFigure);
            app.LowGainLabel.Position = [170 190 200 22];
            app.LowGainLabel.Text = 'Low Band Gain';

            % Create MidGainSlider
            app.MidGainSlider = uislider(app.UIFigure);
            app.MidGainSlider.ValueChangedFcn = createCallbackFcn(app, @MidGainSliderValueChanged, true);
            app.MidGainSlider.Position = [200 140 200 3];
            app.MidGainSlider.Limits = [-20 20]; % Increase range for stronger effect
            app.MidGainSlider.Value = 0; % Set initial value
            app.MidGainSlider.MajorTicks = -20:2:20; % Set major ticks
            app.MidGainSlider.MinorTicks = -20:0.5:20; % Set minor ticks

            % Create MidGainLabel
            app.MidGainLabel = uilabel(app.UIFigure);
            app.MidGainLabel.Position = [170 140 200 22];
            app.MidGainLabel.Text = 'Mid Band Gain';

            % Create HighGainSlider
            app.HighGainSlider = uislider(app.UIFigure);
            app.HighGainSlider.ValueChangedFcn = createCallbackFcn(app, @HighGainSliderValueChanged, true);
            app.HighGainSlider.Position = [200 90 200 3];
            app.HighGainSlider.Limits = [-20 20]; % Increase range for stronger effect
            app.HighGainSlider.Value = 0; % Set initial value
            app.HighGainSlider.MajorTicks = -20:2:20; % Set major ticks
            app.HighGainSlider.MinorTicks = -20:0.5:20; % Set minor ticks

            % Create HighGainLabel
            app.HighGainLabel = uilabel(app.UIFigure);
            app.HighGainLabel.Position = [170 90 200 22];
            app.HighGainLabel.Text = 'High Band Gain';

            % Create WetnessSlider
            app.WetnessSlider = uislider(app.UIFigure);
            app.WetnessSlider.ValueChangedFcn = createCallbackFcn(app, @WetnessSliderValueChanged, true);
            app.WetnessSlider.Position = [200 30 200 3];
            app.WetnessSlider.Limits = [0 100]; % Set minimum and maximum values
            app.WetnessSlider.Value = 50; % Set initial value
            app.WetnessSlider.MajorTicks = 0:10:100; % Set major ticks
            app.WetnessSlider.MinorTicks = 0:1:100; % Set minor ticks

            % Create WetnessLabel
            app.WetnessLabel = uilabel(app.UIFigure);
            app.WetnessLabel.Position = [170 30 200 22];
            app.WetnessLabel.Text = 'Wetness';

            % Create OriginalAxes
            app.OriginalAxes = uiaxes(app.UIFigure);
            title(app.OriginalAxes, 'Original Audio')
            xlabel(app.OriginalAxes, 'Frequency (kHz)')
            ylabel(app.OriginalAxes, 'Power (dB)')
            app.OriginalAxes.Position = [170 250 300 150];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App initialization and construction
    methods (Access = public)

        % Construct app
        function app = AudioFilterApp

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end

% Helper function to design shelf filters
function [b, a] = design_shelf_filter(fs, fc, gain, type)
    A = sqrt(10^(gain/20));
    omega = 2 * pi * fc / fs;
    alpha = sin(omega) / 2 * sqrt((A + 1/A) * (1/0.707 - 1) + 2);
    if strcmp(type, 'low')
        b0 = A * ((A+1) - (A-1)*cos(omega) + 2*sqrt(A)*alpha);
        b1 = 2 * A * ((A-1) - (A+1)*cos(omega));
        b2 = A * ((A+1) - (A-1)*cos(omega) - 2*sqrt(A)*alpha);
        a0 = (A+1) + (A-1)*cos(omega) + 2*sqrt(A)*alpha;
        a1 = -2 * ((A-1) + (A+1)*cos(omega));
        a2 = (A+1) + (A-1)*cos(omega) - 2*sqrt(A)*alpha;
    elseif strcmp(type, 'high')
        b0 = A * ((A+1) + (A-1)*cos(omega) + 2*sqrt(A)*alpha);
        b1 = -2 * A * ((A-1) + (A+1)*cos(omega));
        b2 = A * ((A+1) + (A-1)*cos(omega) - 2*sqrt(A)*alpha);
        a0 = (A+1) - (A-1)*cos(omega) + 2*sqrt(A)*alpha;
        a1 = 2 * ((A-1) - (A+1)*cos(omega));
        a2 = (A+1) - (A-1)*cos(omega) - 2*sqrt(A)*alpha;
    end
    b = [b0 b1 b2] / a0;
    a = [a0 a1 a2] / a0;
end

% Helper function to design peaking filter
function [b, a] = design_peak_filter(fs, fc, gain, Q)
    A = sqrt(10^(gain/20));
    omega = 2 * pi * fc / fs;
    alpha = sin(omega) / (2*Q);
    b0 = 1 + alpha * A;
    b1 = -2 * cos(omega);
    b2 = 1 - alpha * A;
    a0 = 1 + alpha / A;
    a1 = -2 * cos(omega);
    a2 = 1 - alpha / A;
    b = [b0 b1 b2] / a0;
    a = [a0 a1 a2] / a0;
end