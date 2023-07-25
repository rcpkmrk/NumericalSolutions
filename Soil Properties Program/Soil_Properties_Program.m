classdef Soil_Properties_Program < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        MassMeasurementsPanel          matlab.ui.container.Panel
        VoidRatioEditField             matlab.ui.control.NumericEditField
        VoidRatioLabel                 matlab.ui.control.Label
        GsEditField                    matlab.ui.control.NumericEditField
        GsLabel                        matlab.ui.control.Label
        DryMassgEditField              matlab.ui.control.NumericEditField
        DryMassgLabel                  matlab.ui.control.Label
        BulkMassgEditField             matlab.ui.control.NumericEditField
        BulkMassgLabel                 matlab.ui.control.Label
        ResultsPanel                   matlab.ui.container.Panel
        EditField_2                    matlab.ui.control.EditField
        EditField                      matlab.ui.control.EditField
        PlasticityIndexEditField       matlab.ui.control.NumericEditField
        PlasticityIndexEditFieldLabel  matlab.ui.control.Label
        RelativeDensityEditField       matlab.ui.control.NumericEditField
        RelativeDensityEditFieldLabel  matlab.ui.control.Label
        WaterContentEditField          matlab.ui.control.NumericEditField
        WaterContentLabel              matlab.ui.control.Label
        GranularSoilsOnlyPanel         matlab.ui.container.Panel
        e_maxEditField                 matlab.ui.control.NumericEditField
        e_maxLabel                     matlab.ui.control.Label
        e_minEditField                 matlab.ui.control.NumericEditField
        e_minLabel                     matlab.ui.control.Label
        CohesiveSoilsOnlyPanel         matlab.ui.container.Panel
        LiquidLimitEditField           matlab.ui.control.NumericEditField
        LiquidLimitEditFieldLabel      matlab.ui.control.Label
        PlasticLimitEditField          matlab.ui.control.NumericEditField
        PlasticLimitLabel              matlab.ui.control.Label
        ClearAll                       matlab.ui.control.Button
        Calculate                      matlab.ui.control.Button
        SoilTypeSelectionButtonGroup   matlab.ui.container.ButtonGroup
        Cohesive                       matlab.ui.control.RadioButton
        Granular                       matlab.ui.control.RadioButton
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: Calculate
        function CalculatePushed(app, event)
            app.ResultsPanel.Visible = "on"; %results panel appear when calculate is pushed
            app.WaterContentEditField.Value = (app.BulkMassgEditField.Value-app.DryMassgEditField.Value)/...
                (app.DryMassgEditField.Value)*100;                           %water content calculation           
            if (app.Granular.Value == 1)
                app.RelativeDensityEditField.Visible = "on";
                app.RelativeDensityEditField.Value = (app.e_maxEditField.Value-app.VoidRatioEditField.Value)/...
                    (app.e_maxEditField.Value-app.e_minEditField.Value)*100; %relative density calculation 
                app.PlasticityIndexEditField.Visible = "off";
                app.EditField_2.Visible = "on";
                app.EditField_2.Value = "N/A";                               %plasticity index turns N/A
                app.EditField.Visible = "off"; 
            end
            if (app.Cohesive.Value == 1)
                app.PlasticityIndexEditField.Visible = "on";
                app.PlasticityIndexEditField.Value = app.LiquidLimitEditField.Value-...
                    app.PlasticLimitEditField.Value;                         %plasticity index calculation
                app.RelativeDensityEditField.Visible = "off";
                app.EditField.Visible = "on";
                app.EditField.Value = "N/A";                                 %relative density turns N/A
                app.EditField_2.Visible = "off";
            end    
        end

        % Selection changed function: SoilTypeSelectionButtonGroup
        function SoilTypeSelectionButtonGroupSelectionChanged(app, event)
            selectedButton = app.SoilTypeSelectionButtonGroup.SelectedObject;
            if (app.Granular.Value == 1)
                app.GranularSoilsOnlyPanel.Visible = "on"; %if granular chosen, show e_max-e_min panel
                app.CohesiveSoilsOnlyPanel.Visible = "off";
                
            end     
            if (app.Cohesive.Value == 1)
                app.CohesiveSoilsOnlyPanel.Visible = "on"; %if cohesive chosen, show plastic-liquid limit panel
                app.GranularSoilsOnlyPanel.Visible = "off";
                
            end
        end

        % Button pushed function: ClearAll
        function ClearAllPushed(app, event)
            app.ResultsPanel.Visible = "off";
            app.BulkMassgEditField.Value = 0;
            app.DryMassgEditField.Value = 0;
            app.GsEditField.Value = 0;
            app.VoidRatioEditField.Value = 0;
            app.e_maxEditField.Value = 0;
            app.e_minEditField.Value = 0;
            app.PlasticLimitEditField.Value = 0;
            app.LiquidLimitEditField.Value = 0;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 606 349];
            app.UIFigure.Name = 'MATLAB App';

            % Create SoilTypeSelectionButtonGroup
            app.SoilTypeSelectionButtonGroup = uibuttongroup(app.UIFigure);
            app.SoilTypeSelectionButtonGroup.SelectionChangedFcn = createCallbackFcn(app, ...
                @SoilTypeSelectionButtonGroupSelectionChanged, true);
            app.SoilTypeSelectionButtonGroup.Title = 'Soil Type Selection';
            app.SoilTypeSelectionButtonGroup.Position = [30 264 123 72];

            % Create Granular
            app.Granular = uiradiobutton(app.SoilTypeSelectionButtonGroup);
            app.Granular.Text = 'Granular';
            app.Granular.Position = [11 26 69 22];
            app.Granular.Value = true;

            % Create Cohesive
            app.Cohesive = uiradiobutton(app.SoilTypeSelectionButtonGroup);
            app.Cohesive.Text = 'Cohesive';
            app.Cohesive.Position = [11 4 72 22];

            % Create Calculate
            app.Calculate = uibutton(app.UIFigure, 'push');
            app.Calculate.ButtonPushedFcn = createCallbackFcn(app, @CalculatePushed, true);
            app.Calculate.BackgroundColor = [0.302 0.7451 0.9333];
            app.Calculate.Position = [306 184 100 22];
            app.Calculate.Text = 'Calculate';

            % Create ClearAll
            app.ClearAll = uibutton(app.UIFigure, 'push');
            app.ClearAll.ButtonPushedFcn = createCallbackFcn(app, @ClearAllPushed, true);
            app.ClearAll.BackgroundColor = [0.302 0.7451 0.9333];
            app.ClearAll.Position = [464 184 100 22];
            app.ClearAll.Text = 'Clear All';

            % Create CohesiveSoilsOnlyPanel
            app.CohesiveSoilsOnlyPanel = uipanel(app.UIFigure);
            app.CohesiveSoilsOnlyPanel.Title = 'Cohesive Soils Only';
            app.CohesiveSoilsOnlyPanel.Visible = 'off';
            app.CohesiveSoilsOnlyPanel.Position = [30 26 179 93];

            % Create PlasticLimitLabel
            app.PlasticLimitLabel = uilabel(app.CohesiveSoilsOnlyPanel);
            app.PlasticLimitLabel.HorizontalAlignment = 'right';
            app.PlasticLimitLabel.Position = [-8 47 175 22];
            app.PlasticLimitLabel.Text = 'Plastic Limit:                        %';

            % Create PlasticLimitEditField
            app.PlasticLimitEditField = uieditfield(app.CohesiveSoilsOnlyPanel, 'numeric');
            app.PlasticLimitEditField.HorizontalAlignment = 'center';
            app.PlasticLimitEditField.Tooltip = {''};
            app.PlasticLimitEditField.Position = [81 47 65 22];
            app.PlasticLimitEditField.Value = 0;

            % Create LiquidLimitEditFieldLabel
            app.LiquidLimitEditFieldLabel = uilabel(app.CohesiveSoilsOnlyPanel);
            app.LiquidLimitEditFieldLabel.HorizontalAlignment = 'right';
            app.LiquidLimitEditFieldLabel.Position = [-7 15 175 22];
            app.LiquidLimitEditFieldLabel.Text = 'Liquid Limit:                        %';

            % Create LiquidLimitEditField
            app.LiquidLimitEditField = uieditfield(app.CohesiveSoilsOnlyPanel, 'numeric');
            app.LiquidLimitEditField.HorizontalAlignment = 'center';
            app.LiquidLimitEditField.Position = [82 15 65 22];
            app.LiquidLimitEditField.Value = 0;

            % Create GranularSoilsOnlyPanel
            app.GranularSoilsOnlyPanel = uipanel(app.UIFigure);
            app.GranularSoilsOnlyPanel.Title = 'Granular Soils Only';
            app.GranularSoilsOnlyPanel.Position = [30 141 143 97];

            % Create e_minLabel
            app.e_minLabel = uilabel(app.GranularSoilsOnlyPanel);
            app.e_minLabel.HorizontalAlignment = 'center';
            app.e_minLabel.Position = [11 51 43 22];
            app.e_minLabel.Text = 'e_min:';

            % Create e_minEditField
            app.e_minEditField = uieditfield(app.GranularSoilsOnlyPanel, 'numeric');
            app.e_minEditField.HorizontalAlignment = 'center';
            app.e_minEditField.Position = [60 51 67 22];
            app.e_minEditField.Value = 0;

            % Create e_maxLabel
            app.e_maxLabel = uilabel(app.GranularSoilsOnlyPanel);
            app.e_maxLabel.HorizontalAlignment = 'right';
            app.e_maxLabel.Position = [7 21 48 22];
            app.e_maxLabel.Text = 'e_max: ';

            % Create e_maxEditField
            app.e_maxEditField = uieditfield(app.GranularSoilsOnlyPanel, 'numeric');
            app.e_maxEditField.HorizontalAlignment = 'center';
            app.e_maxEditField.Position = [60 21 68 22];
            app.e_maxEditField.Value = 0;

            % Create ResultsPanel
            app.ResultsPanel = uipanel(app.UIFigure);
            app.ResultsPanel.TitlePosition = 'centertop';
            app.ResultsPanel.Title = 'Results';
            app.ResultsPanel.Visible = 'off';
            app.ResultsPanel.FontWeight = 'bold';
            app.ResultsPanel.Position = [306 26 257 141];

            % Create WaterContentLabel
            app.WaterContentLabel = uilabel(app.ResultsPanel);
            app.WaterContentLabel.HorizontalAlignment = 'right';
            app.WaterContentLabel.FontWeight = 'bold';
            app.WaterContentLabel.Position = [11 94 91 22];
            app.WaterContentLabel.Text = 'Water Content:';

            % Create WaterContentEditField
            app.WaterContentEditField = uieditfield(app.ResultsPanel, 'numeric');
            app.WaterContentEditField.Editable = 'off';
            app.WaterContentEditField.HorizontalAlignment = 'center';
            app.WaterContentEditField.FontWeight = 'bold';
            app.WaterContentEditField.BackgroundColor = [0.9412 0.9412 0.9412];
            app.WaterContentEditField.Position = [117 94 100 22];

            % Create RelativeDensityEditFieldLabel
            app.RelativeDensityEditFieldLabel = uilabel(app.ResultsPanel);
            app.RelativeDensityEditFieldLabel.HorizontalAlignment = 'right';
            app.RelativeDensityEditFieldLabel.FontWeight = 'bold';
            app.RelativeDensityEditFieldLabel.Position = [0 55 102 22];
            app.RelativeDensityEditFieldLabel.Text = 'Relative Density:';

            % Create RelativeDensityEditField
            app.RelativeDensityEditField = uieditfield(app.ResultsPanel, 'numeric');
            app.RelativeDensityEditField.Editable = 'off';
            app.RelativeDensityEditField.HorizontalAlignment = 'center';
            app.RelativeDensityEditField.FontWeight = 'bold';
            app.RelativeDensityEditField.BackgroundColor = [0.9412 0.9412 0.9412];
            app.RelativeDensityEditField.Position = [117 55 100 22];

            % Create PlasticityIndexEditFieldLabel
            app.PlasticityIndexEditFieldLabel = uilabel(app.ResultsPanel);
            app.PlasticityIndexEditFieldLabel.HorizontalAlignment = 'right';
            app.PlasticityIndexEditFieldLabel.FontWeight = 'bold';
            app.PlasticityIndexEditFieldLabel.Position = [5 17 97 22];
            app.PlasticityIndexEditFieldLabel.Text = 'Plasticity Index:';

            % Create PlasticityIndexEditField
            app.PlasticityIndexEditField = uieditfield(app.ResultsPanel, 'numeric');
            app.PlasticityIndexEditField.Editable = 'off';
            app.PlasticityIndexEditField.HorizontalAlignment = 'center';
            app.PlasticityIndexEditField.FontWeight = 'bold';
            app.PlasticityIndexEditField.BackgroundColor = [0.9412 0.9412 0.9412];
            app.PlasticityIndexEditField.Position = [117 17 100 22];

            % Create EditField
            app.EditField = uieditfield(app.ResultsPanel, 'text');
            app.EditField.Editable = 'off';
            app.EditField.HorizontalAlignment = 'center';
            app.EditField.FontWeight = 'bold';
            app.EditField.BackgroundColor = [0.9412 0.9412 0.9412];
            app.EditField.Position = [117 55 100 22];

            % Create EditField_2
            app.EditField_2 = uieditfield(app.ResultsPanel, 'text');
            app.EditField_2.Editable = 'off';
            app.EditField_2.HorizontalAlignment = 'center';
            app.EditField_2.FontWeight = 'bold';
            app.EditField_2.BackgroundColor = [0.9412 0.9412 0.9412];
            app.EditField_2.Position = [117 17 100 22];

            % Create MassMeasurementsPanel
            app.MassMeasurementsPanel = uipanel(app.UIFigure);
            app.MassMeasurementsPanel.Title = 'Mass Measurements';
            app.MassMeasurementsPanel.Position = [249 229 329 107];

            % Create BulkMassgLabel
            app.BulkMassgLabel = uilabel(app.MassMeasurementsPanel);
            app.BulkMassgLabel.HorizontalAlignment = 'right';
            app.BulkMassgLabel.Position = [12 55 77 22];
            app.BulkMassgLabel.Text = 'Bulk Mass (g):';

            % Create BulkMassgEditField
            app.BulkMassgEditField = uieditfield(app.MassMeasurementsPanel, 'numeric');
            app.BulkMassgEditField.HorizontalAlignment = 'center';
            app.BulkMassgEditField.Position = [99 55 59 26];
            app.BulkMassgEditField.Value = 0;

            % Create DryMassgLabel
            app.DryMassgLabel = uilabel(app.MassMeasurementsPanel);
            app.DryMassgLabel.HorizontalAlignment = 'right';
            app.DryMassgLabel.Position = [12 11 77 22];
            app.DryMassgLabel.Text = 'Dry Mass (g):';

            % Create DryMassgEditField
            app.DryMassgEditField = uieditfield(app.MassMeasurementsPanel, 'numeric');
            app.DryMassgEditField.HorizontalAlignment = 'center';
            app.DryMassgEditField.Position = [99 11 59 26];
            app.DryMassgEditField.Value = 0;

            % Create GsLabel
            app.GsLabel = uilabel(app.MassMeasurementsPanel);
            app.GsLabel.HorizontalAlignment = 'right';
            app.GsLabel.Position = [213 53 34 22];
            app.GsLabel.Text = 'Gs:';

            % Create GsEditField
            app.GsEditField = uieditfield(app.MassMeasurementsPanel, 'numeric');
            app.GsEditField.HorizontalAlignment = 'center';
            app.GsEditField.Position = [256 53 59 26];
            app.GsEditField.Value = 0;

            % Create VoidRatioLabel
            app.VoidRatioLabel = uilabel(app.MassMeasurementsPanel);
            app.VoidRatioLabel.HorizontalAlignment = 'right';
            app.VoidRatioLabel.Position = [185 11 63 22];
            app.VoidRatioLabel.Text = 'Void Ratio:';

            % Create VoidRatioEditField
            app.VoidRatioEditField = uieditfield(app.MassMeasurementsPanel, 'numeric');
            app.VoidRatioEditField.HorizontalAlignment = 'center';
            app.VoidRatioEditField.Position = [257 11 59 26];
            app.VoidRatioEditField.Value = 0;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Soil_Properties_Program

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