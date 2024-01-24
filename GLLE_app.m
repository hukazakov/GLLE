classdef GLLE_app < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        cNLField                matlab.ui.control.EditField
        Label_6                 matlab.ui.control.Label
        cDField                 matlab.ui.control.EditField
        Label_5                 matlab.ui.control.Label
        ParameterButtonGroup    matlab.ui.container.ButtonGroup
        GVDLabel                matlab.ui.control.Label
        PumpLabel               matlab.ui.control.Label
        KerrLabel               matlab.ui.control.Label
        LEFLabel                matlab.ui.control.Label
        xiField                 matlab.ui.control.EditField
        Label_4                 matlab.ui.control.Label
        gammaField              matlab.ui.control.EditField
        Label_3                 matlab.ui.control.Label
        betaField               matlab.ui.control.EditField
        Label_2                 matlab.ui.control.Label
        alphaField              matlab.ui.control.EditField
        Label                   matlab.ui.control.Label
        ModelButtonGroup        matlab.ui.container.ButtonGroup
        LLEButton               matlab.ui.control.ToggleButton
        CGLEButton              matlab.ui.control.ToggleButton
        FCGLEButton             matlab.ui.control.ToggleButton
        RunButton               matlab.ui.control.StateButton
        SavestateButton         matlab.ui.control.Button
        ThetaField              matlab.ui.control.NumericEditField
        ThetaSlider             matlab.ui.control.Slider
        SliderLabel             matlab.ui.control.Label
        YField                  matlab.ui.control.NumericEditField
        YSlider                 matlab.ui.control.Slider
        YSliderLabel            matlab.ui.control.Label
        SpectrumUIAxes          matlab.ui.control.UIAxes
        CurrentRoundtripUIAxes  matlab.ui.control.UIAxes
        MeanFieldUIAxes         matlab.ui.control.UIAxes
        SpatiotemporalUIAxes    matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        s = FCGLE_simulator();% Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.s = app.s.initialize_solver();
            app.s.UI_axis_cavity_intensity = app.CurrentRoundtripUIAxes;
            app.s.UI_axis_mean_intensity = app.MeanFieldUIAxes;
            app.s.UI_axis_spacetime = app.SpatiotemporalUIAxes;
            app.s.UI_axis_spectrum = app.SpectrumUIAxes;
            
            app.ThetaField.Value = app.s.Theta;
            app.ThetaSlider.Value = app.s.Theta;
            app.YField.Value = app.s.Y;
            app.YSlider.Value = app.s.Y;
        end

        % Callback function
        function RunButtonPushed(app, event)

        end

        % Callback function
        function ParameterDropDownValueChanged(app, event)
            value = app.ParameterDropDown.Value;
            
            
            switch value
                case 'α'
                    app.ParameterField.Value = app.s.alpha;
                case 'β'
                    app.ParameterField.Value = app.s.beta;
                case 'Γ'
                    app.ParameterField.Value = app.s.gamma;
            end
        end

        % Callback function
        function ThetaFieldValueChanged(app, event)
            value = app.ThetaField.Value;
            app.ThetaSlider.Value = value;
            
            app.s.run_simulation = false;
            app.s.Theta = value;
            app.s.run_simulation = true;
            
        end

        % Value changing function: YSlider
        function YSliderValueChanging(app, event)
            changingValue = event.Value;
            app.YField.Value = changingValue;
            
            app.s.run_simulation = false;
            app.s.Y = changingValue;
            app.s.run_simulation = true;
        end

        % Value changed function: YField
        function YFieldValueChanged(app, event)
            value = app.YField.Value;
            app.YSlider.Value = value;
            
            app.s.run_simulation = false;
            app.s.Y = value;
            app.s.run_simulation = true;
        end

        % Value changing function: ThetaSlider
        function ThetaSliderValueChanging(app, event)
            changingValue = event.Value;
            app.ThetaField.Value = changingValue;
            
            app.s.run_simulation = false;
            app.s.Theta = changingValue;
            app.s.run_simulation = true;
        end

        % Value changed function: ThetaField
        function ThetaFieldValueChanged2(app, event)
            value = app.ThetaField.Value;
            app.ThetaSlider.Value = value;
        end

        % Callback function
        function ParameterFieldValueChanged(app, event)
           
        end

        % Value changed function: RunButton
        function RunButtonValueChanged(app, event)
            while 1
                value = app.RunButton.Value;
                switch value
                    case 1
                        app.RunButton.Text = "Pause";
                        app.RunButton.BackgroundColor = [1.00,0.71,0.57];
                        
                        app.s.Y = app.YField.Value;
                        app.s.Theta = app.ThetaField.Value;
                        app.s = app.s.update_params();

                        app.s.run_simulation = true;
                        app.s = app.s.advance_simulation;
                        pause(0.05)
                    case 0
                        app.s.run_simulation = false;
                        app.RunButton.Text = "Run";
                        app.RunButton.BackgroundColor = [0.31,0.78,0.47];
                        break
                end
      
            end
            
        end

        % Selection changed function: ModelButtonGroup
        function ModelButtonGroupSelectionChanged(app, event)
            selectedButton = app.ModelButtonGroup.SelectedObject;
            switch selectedButton.Text
                case "FCGLE"
                    app.alphaField.Value = 2;
                    app.betaField.Value = 0;
                    app.gammaField.Value = 1;
                    app.xiField.Value = 1;

                    app.YField.Value = 6;
                    app.YSlider.Value = 6;

                    app.cNLField.Value = 0;
                    app.cNLField.Enable = false;

                    app.cDField.Value = 0;
                    app.cDField.Enable = false;

                    app.s.u = zeros(size(app.s.x));
                case "CGLE"
                    app.alphaField.Value = 2;
                    app.betaField.Value = 0;
                    app.gammaField.Value = 1;
                    app.xiField.Value = 1;

                    app.YField.Value = 0;
                    app.YSlider.Value = 0;

                    app.cNLField.Value = 1;
                    app.cNLField.Enable = true;
                    app.cDField.Value = 1;
                    app.cDField.Enable = true;
                    app.cDField.FontColor = [0 0 0];

                    app.s.u = zeros(size(app.s.x));
                case "LLE"
                    app.alphaField.Value = 0;
%                     app.betaField.Value = 8;
                    app.betaField.Value = "-1+1i";
                    app.gammaField.Value = -1;
%                     app.xiField.Value = 4;
                    app.xiField.Value = "1+1i";

                    app.YField.Value = 6;
                    app.YSlider.Value = 6;  

                    app.cNLField.Value = 0;
                    app.cNLField.Enable = false;

                    app.cDField.Value = 0;
                    app.cDField.Enable = false;

                    app.betaField.Enable = true;
                    app.xiField.Enable = true;
                    
                    % reset waveform when switching to LLE
                    app.s.u = zeros(size(app.s.x));
            end

            app.s.alpha = app.alphaField.Value;
            app.s.beta = app.betaField.Value;
            app.s.gamma = app.gammaField.Value;
            app.s.xi = app.xiField.Value;
            app.s.Y = app.YField.Value;

            app.s = app.s.update_params();
        end

        % Value changed function: alphaField
        function alphaFieldValueChanged(app, event)
            value = app.alphaField.Value;
            app.s.alpha = value;
            app.s = app.s.update_params();
        end

        % Value changed function: betaField
        function betaFieldValueChanged(app, event)
            value = app.betaField.Value;
            app.s.beta = value;
            app.s = app.s.update_params();
            
        end

        % Value changed function: gammaField
        function gammaFieldValueChanged(app, event)
            value = app.gammaField.Value;
            app.s.gamma = value;
            app.s = app.s.update_params();
        end

        % Value changed function: xiField
        function xiFieldValueChanged(app, event)
            value = app.xiField.Value;
            app.s.xi = value;
            app.s = app.s.update_params();
        end

        % Button pushed function: SavestateButton
        function SavestateButtonPushed(app, event)
%             vars = {'app.s.u_spatiotemporal', 'app.s.u_spectrum'};
%             params = {app.s.alpha, app.s.beta, app.s.gamma, app.s.;
% 
%             uisave(vars)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1500 750];
            app.UIFigure.Name = 'MATLAB App';

            % Create SpatiotemporalUIAxes
            app.SpatiotemporalUIAxes = uiaxes(app.UIFigure);
            ylabel(app.SpatiotemporalUIAxes, 'Slow time')
            zlabel(app.SpatiotemporalUIAxes, 'Z')
            app.SpatiotemporalUIAxes.XLim = [-50 50];
            app.SpatiotemporalUIAxes.YLim = [0 300];
            app.SpatiotemporalUIAxes.XColor = 'none';
            app.SpatiotemporalUIAxes.XTick = [];
            app.SpatiotemporalUIAxes.YAxisLocation = 'right';
            app.SpatiotemporalUIAxes.YTick = [0 300];
            app.SpatiotemporalUIAxes.YTickLabel = {'0'; '300'};
            app.SpatiotemporalUIAxes.TickDir = 'out';
            app.SpatiotemporalUIAxes.Position = [942 257 540 467];

            % Create MeanFieldUIAxes
            app.MeanFieldUIAxes = uiaxes(app.UIFigure);
            xlabel(app.MeanFieldUIAxes, 'Slow time')
            ylabel(app.MeanFieldUIAxes, 'Mean intensity')
            zlabel(app.MeanFieldUIAxes, 'Z')
            app.MeanFieldUIAxes.YLim = [0 3];
            app.MeanFieldUIAxes.XTickLabel = '';
            app.MeanFieldUIAxes.YTick = [0 5];
            app.MeanFieldUIAxes.YTickLabel = {'0'; '5'};
            app.MeanFieldUIAxes.TickDir = 'out';
            app.MeanFieldUIAxes.Position = [396 25 527 228];

            % Create CurrentRoundtripUIAxes
            app.CurrentRoundtripUIAxes = uiaxes(app.UIFigure);
            xlabel(app.CurrentRoundtripUIAxes, 'Cavity position')
            ylabel(app.CurrentRoundtripUIAxes, 'Intensity')
            zlabel(app.CurrentRoundtripUIAxes, 'Z')
            app.CurrentRoundtripUIAxes.XLim = [-50 50];
            app.CurrentRoundtripUIAxes.YLim = [0 5];
            app.CurrentRoundtripUIAxes.XTickLabel = '';
            app.CurrentRoundtripUIAxes.YAxisLocation = 'right';
            app.CurrentRoundtripUIAxes.YTick = [0 5];
            app.CurrentRoundtripUIAxes.TickDir = 'out';
            app.CurrentRoundtripUIAxes.Position = [942 25 527 227];

            % Create SpectrumUIAxes
            app.SpectrumUIAxes = uiaxes(app.UIFigure);
            xlabel(app.SpectrumUIAxes, 'Frequency (FSR)')
            ylabel(app.SpectrumUIAxes, 'PSD (dB)')
            zlabel(app.SpectrumUIAxes, 'Z')
            app.SpectrumUIAxes.TickDir = 'out';
            app.SpectrumUIAxes.Position = [386 258 537 465];

            % Create YSliderLabel
            app.YSliderLabel = uilabel(app.UIFigure);
            app.YSliderLabel.HorizontalAlignment = 'center';
            app.YSliderLabel.FontSize = 24;
            app.YSliderLabel.Position = [27 517 25 31];
            app.YSliderLabel.Text = 'Y';

            % Create YSlider
            app.YSlider = uislider(app.UIFigure);
            app.YSlider.Limits = [0 20];
            app.YSlider.MajorTicks = [0 5 10 15 20];
            app.YSlider.ValueChangingFcn = createCallbackFcn(app, @YSliderValueChanging, true);
            app.YSlider.Position = [73 541 150 3];

            % Create YField
            app.YField = uieditfield(app.UIFigure, 'numeric');
            app.YField.ValueChangedFcn = createCallbackFcn(app, @YFieldValueChanged, true);
            app.YField.FontSize = 24;
            app.YField.Position = [253 523 100 31];

            % Create SliderLabel
            app.SliderLabel = uilabel(app.UIFigure);
            app.SliderLabel.HorizontalAlignment = 'center';
            app.SliderLabel.FontSize = 24;
            app.SliderLabel.Position = [27 581 25 31];
            app.SliderLabel.Text = 'θ';

            % Create ThetaSlider
            app.ThetaSlider = uislider(app.UIFigure);
            app.ThetaSlider.Limits = [-10 10];
            app.ThetaSlider.MajorTicks = [-10 -5 0 5 10];
            app.ThetaSlider.ValueChangingFcn = createCallbackFcn(app, @ThetaSliderValueChanging, true);
            app.ThetaSlider.Position = [73 605 150 3];

            % Create ThetaField
            app.ThetaField = uieditfield(app.UIFigure, 'numeric');
            app.ThetaField.ValueChangedFcn = createCallbackFcn(app, @ThetaFieldValueChanged2, true);
            app.ThetaField.FontSize = 24;
            app.ThetaField.Position = [253 587 100 31];

            % Create SavestateButton
            app.SavestateButton = uibutton(app.UIFigure, 'push');
            app.SavestateButton.ButtonPushedFcn = createCallbackFcn(app, @SavestateButtonPushed, true);
            app.SavestateButton.FontSize = 24;
            app.SavestateButton.Position = [202 649 151 48];
            app.SavestateButton.Text = 'Save state';

            % Create RunButton
            app.RunButton = uibutton(app.UIFigure, 'state');
            app.RunButton.ValueChangedFcn = createCallbackFcn(app, @RunButtonValueChanged, true);
            app.RunButton.Text = 'Run';
            app.RunButton.BackgroundColor = [0.3098 0.7804 0.4706];
            app.RunButton.FontSize = 24;
            app.RunButton.Position = [73 649 117 48];

            % Create ModelButtonGroup
            app.ModelButtonGroup = uibuttongroup(app.UIFigure);
            app.ModelButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ModelButtonGroupSelectionChanged, true);
            app.ModelButtonGroup.BorderType = 'none';
            app.ModelButtonGroup.TitlePosition = 'centertop';
            app.ModelButtonGroup.Title = 'Model';
            app.ModelButtonGroup.FontSize = 18;
            app.ModelButtonGroup.Position = [73 249 123 216];

            % Create FCGLEButton
            app.FCGLEButton = uitogglebutton(app.ModelButtonGroup);
            app.FCGLEButton.Text = 'FCGLE';
            app.FCGLEButton.FontSize = 18;
            app.FCGLEButton.Position = [6 143 113 38];
            app.FCGLEButton.Value = true;

            % Create CGLEButton
            app.CGLEButton = uitogglebutton(app.ModelButtonGroup);
            app.CGLEButton.Text = 'CGLE';
            app.CGLEButton.FontSize = 18;
            app.CGLEButton.Position = [6 97 113 38];

            % Create LLEButton
            app.LLEButton = uitogglebutton(app.ModelButtonGroup);
            app.LLEButton.Text = 'LLE';
            app.LLEButton.FontSize = 18;
            app.LLEButton.Position = [6 51 113 38];

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.HorizontalAlignment = 'center';
            app.Label.FontSize = 18;
            app.Label.Position = [223 398 25 23];
            app.Label.Text = 'α';

            % Create alphaField
            app.alphaField = uieditfield(app.UIFigure, 'numeric');
            app.alphaField.ValueChangedFcn = createCallbackFcn(app, @alphaFieldValueChanged, true);
            app.alphaField.FontSize = 18;
            app.alphaField.Position = [259 391 39 37];
            app.alphaField.Value = 2;

            % Create Label_2
            app.Label_2 = uilabel(app.UIFigure);
            app.Label_2.HorizontalAlignment = 'center';
            app.Label_2.FontSize = 18;
            app.Label_2.Enable = 'off';
            app.Label_2.Position = [223 349 25 23];
            app.Label_2.Text = 'β';

            % Create betaField
            app.betaField = uieditfield(app.UIFigure, 'numeric');
            app.betaField.ValueChangedFcn = createCallbackFcn(app, @betaFieldValueChanged, true);
            app.betaField.FontSize = 18;
            app.betaField.Enable = 'off';
            app.betaField.Position = [259 342 39 37];

            % Create Label_3
            app.Label_3 = uilabel(app.UIFigure);
            app.Label_3.HorizontalAlignment = 'center';
            app.Label_3.FontSize = 18;
            app.Label_3.Enable = 'off';
            app.Label_3.Position = [223 300 25 23];
            app.Label_3.Text = 'γ';

            % Create gammaField
            app.gammaField = uieditfield(app.UIFigure, 'numeric');
            app.gammaField.ValueChangedFcn = createCallbackFcn(app, @gammaFieldValueChanged, true);
            app.gammaField.FontSize = 18;
            app.gammaField.Enable = 'off';
            app.gammaField.Position = [259 293 39 37];
            app.gammaField.Value = 1;

            % Create Label_4
            app.Label_4 = uilabel(app.UIFigure);
            app.Label_4.HorizontalAlignment = 'center';
            app.Label_4.FontSize = 18;
            app.Label_4.Position = [223 251 25 23];
            app.Label_4.Text = 'ζ';

            % Create xiField
            app.xiField = uieditfield(app.UIFigure, 'numeric');
            app.xiField.ValueChangedFcn = createCallbackFcn(app, @xiFieldValueChanged, true);
            app.xiField.FontSize = 18;
            app.xiField.Position = [259 244 39 37];
            app.xiField.Value = 1;

            % Create LEFLabel
            app.LEFLabel = uilabel(app.UIFigure);
            app.LEFLabel.FontSize = 18;
            app.LEFLabel.Position = [317 398 37 23];
            app.LEFLabel.Text = 'LEF';

            % Create KerrLabel
            app.KerrLabel = uilabel(app.UIFigure);
            app.KerrLabel.FontSize = 18;
            app.KerrLabel.Position = [317 350 39 23];
            app.KerrLabel.Text = 'Kerr';

            % Create PumpLabel
            app.PumpLabel = uilabel(app.UIFigure);
            app.PumpLabel.FontSize = 18;
            app.PumpLabel.Position = [317 303 53 23];
            app.PumpLabel.Text = 'Pump';

            % Create GVDLabel
            app.GVDLabel = uilabel(app.UIFigure);
            app.GVDLabel.FontSize = 18;
            app.GVDLabel.Position = [317 253 43 23];
            app.GVDLabel.Text = 'GVD';

            % Create ParameterButtonGroup
            app.ParameterButtonGroup = uibuttongroup(app.UIFigure);
            app.ParameterButtonGroup.BorderType = 'none';
            app.ParameterButtonGroup.TitlePosition = 'centertop';
            app.ParameterButtonGroup.Title = 'Parameter';
            app.ParameterButtonGroup.FontSize = 18;
            app.ParameterButtonGroup.Position = [233 435 123 30];

            % Create Label_5
            app.Label_5 = uilabel(app.UIFigure);
            app.Label_5.HorizontalAlignment = 'center';
            app.Label_5.FontSize = 18;
            app.Label_5.Enable = 'off';
            app.Label_5.Position = [222 203 28 23];
            app.Label_5.Text = 'cD';

            % Create cDField
            app.cDField = uieditfield(app.UIFigure, 'numeric');
            app.cDField.FontSize = 18;
            app.cDField.Enable = 'off';
            app.cDField.Position = [259 196 39 37];
            app.cDField.Value = 1;

            % Create Label_6
            app.Label_6 = uilabel(app.UIFigure);
            app.Label_6.HorizontalAlignment = 'center';
            app.Label_6.FontSize = 18;
            app.Label_6.Enable = 'off';
            app.Label_6.Position = [217 155 38 23];
            app.Label_6.Text = 'cNL';

            % Create cNLField
            app.cNLField = uieditfield(app.UIFigure, 'numeric');
            app.cNLField.FontSize = 18;
            app.cNLField.Enable = 'off';
            app.cNLField.Position = [259 148 39 37];
            app.cNLField.Value = 1;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = GLLE_app

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

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