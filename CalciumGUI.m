classdef CalciumGUI < handle
    
    properties
        Name              % name of gui figure
        hFigure           % figure handle
        nPlot             % Number plots
        hSliderPanel      % panel handle containing ui elements
        hButtonPanel      % panel handle containing buttons
        hAxisPanel        % panel handle containing axes
        hButton           % vector of all button handles
        nButton           % number of buttons currently active
        hAxis             % axis handles
        hSlider           % slider handles
        SliderMaxList     % Array containing max slider values
        
        nParam            % number of parameters
        InitialParam      % initial parameter set
        ParamNameList     % list of all parameter names
        ParamValueList    % list of all parameter values
        ParamUnitList     % list of parameter units
        Parameters        % intial slider values
        
        Time              % time vector
        StimOnset         % start of IP3 stimulation
        StateVar          % current state variables
        NonStateVar       % current non-state variables
    end
    
    %% GENERAL
    methods
        
        % Constructor
        function obj = CalciumGUI(P, Name, sliderMaxList, unitList)
            obj.SliderMaxList = sliderMaxList;
            obj.ParamUnitList = unitList;
            obj.Name = Name;
            obj.nPlot = 4;
            obj.hFigure = figure('Visible','on','Name',obj.Name);
            obj.hFigure.Units = 'normalized';
            obj.hFigure.Position = [0.1, 0.1, 0.8, 0.8];
            
            obj.InitialParam = P;
            obj.Parameters = P;
            obj.ParamNameList = fieldnames(obj.Parameters);
            obj.ParamValueList = zeros(obj.nParam,1);
            obj.nParam = numel(obj.ParamNameList);
            
            for i = 1:obj.nParam
                obj.ParamValueList(i) = ...
                    obj.Parameters.(obj.ParamNameList{i});
            end
            
            createSliderPanel(obj)
            createButtonPanel(obj)
            createAxisPanel(obj)
            createAxes(obj)
            setsliders(obj)
            setsliderposition(obj)
            obj.nButton = 0;
            
            addButton(obj, 'RESET', @obj.onReset)
            addButton(obj, 'EXPORT', @obj.onExport)
            addButton(obj, 'IMPORT', @obj.onImport)
            addButton(obj, 'NON-STATE', @obj.onShowNonState)
            addButton(obj, 'BIFURCATION', @obj.onBifDiag)
            addButton(obj, 'NORMALIZED', @obj.onNormalized)
            
            setbuttonposition(obj)
            setinitialaxis(obj)
            updateaxes(obj)
            
        end
        
        % Panel containing parameter sliders
        function createSliderPanel(obj)
            obj.hSliderPanel.Controls = uipanel(obj.hFigure,...
                'Title','Parameters',...
                'FontSize',10,'FontWeight','bold',...
                'Position',[0.6 0.3 0.4 0.7],...
                'BackgroundColor','w',...
                'Units','normalized');
        end
        
        % Panel containing buttons
        function createButtonPanel(obj)
            obj.hButtonPanel.Controls = uipanel(obj.hFigure,...
                'Title','Buttons',...
                'FontSize',10,'FontWeight','bold',...
                'Position',[0.6 0 0.4 0.3],...
                'BackgroundColor','w');
        end
        
        % Panel containing axes
        function createAxisPanel(obj)
            obj.hAxisPanel = uipanel(obj.hFigure,...
                'Position', [0 0 0.6 1],...
                'BackgroundColor','w');
        end
        
        % Create the axes
        function createAxes(obj)
            nrow = 2;
            ncol = 2;
            parent = obj.hAxisPanel;
            % ensure spacing between each axis
            spacing = 0.1;
            width = (1 - (ncol+1)*spacing)/ncol;
            height = (1 - (nrow+1)*spacing)/nrow;
            count = 0;
            for i = 1:nrow
                vert = i*spacing + (i - 1)*height;
                for j = 1:ncol
                    horiz = j*spacing + (j - 1)*width;
                    count = count + 1;
                    pos = [horiz, vert, width, height];
                    obj.hAxis.h(count) = axes(parent, 'Position', pos);
                end
            end
        end
        
        % Define sliders
        function setsliders(obj)
            % loop over all parameters and create sliders
            for i = 1:obj.nParam
                paramName = obj.ParamNameList{i};
                paramValue = obj.ParamValueList(i);
                sMin = 0;
                sMax = obj.SliderMaxList(i);
                obj.hSlider.s(i,1) = ...
                    uicontrol(...
                    'Style','slider',...
                    'String',paramName,...
                    'Min',sMin,'Max',sMax,'Value',paramValue,...
                    'SliderStep',[0.001 0.1],...
                    'Units','normalized',...
                    'Callback',@obj.updateaxes);
                obj.hSlider.t(i,1) = ...
                    uicontrol(...
                    'Style','text',...
                    'Units','normalized');
            end
        end
        
        % Set all slider positions
        function setsliderposition(obj)
            % panel for column 1 of parameter sliders
            col1 = uipanel('Parent',obj.hSliderPanel.Controls,...
                'Position', [0 0 0.5 1],...
                'Units','normalized');
            % panel for column 2 of parameter sliders
            col2 = uipanel('Parent',obj.hSliderPanel.Controls,...
                'Position', [0.5 0 0.5 1],...
                'Units','normalized');
            % parameters go in column 1
            ncol1 = ceil(obj.nParam/2);
            % remaining parameters go in column 2
            ncol2 = obj.nParam - ncol1;
            % loop through column 1 parameters and set positions
            count = 0;
            for i = 1:ncol1
                count = count + 1;
                swidth = 0.4;
                twidth = 1 - swidth;
                height = 1/ncol1;
                spos = [twidth 1-height*i swidth height];
                tpos = [0 1-height*i twidth height];
                obj.hSlider.s(count,1).Position = spos;
                obj.hSlider.t(count,1).Position = tpos;
                obj.hSlider.s(count,1).Parent = col1;
                obj.hSlider.t(count,1).Parent = col1;
                obj.hSlider.t(count,1).HorizontalAlignment = 'left';
                obj.hSlider.t(count,1).FontSize = 4;
            end
            % loop through column 2 parameters and set positions
            for i = 1:ncol2
                count = count + 1;
                swidth = 0.4;
                twidth = 1 - swidth;
                height = 1/ncol2;
                spos = [twidth 1-height*i swidth height];
                tpos = [0 1-height*i twidth height];
                obj.hSlider.s(count,1).Position = spos;
                obj.hSlider.t(count,1).Position = tpos;
                obj.hSlider.s(count,1).Parent = col2;
                obj.hSlider.t(count,1).Parent = col2;
                obj.hSlider.t(count,1).HorizontalAlignment = 'left';
                obj.hSlider.t(count,1).FontSize = 4;
            end
        end
        
        % Update the text showing current slider values
        function updatetext(obj)
            for i = 1:obj.nParam
                name = obj.ParamNameList{i};
                val = obj.hSlider.s(i).Value;
                unit = obj.ParamUnitList{i};
                nsigdig = 3;
                if isempty(unit) % for unit-less parameters
                    sliderText = [name ' ' unit ': ' num2str(round(val,nsigdig,'significant'))];
                else
                    sliderText = [name ' [' unit ']: ' num2str(round(val,nsigdig,'significant'))];
                end
                obj.hSlider.t(i).String = sliderText;
                obj.hSlider.t(i).FontSize = 10;
                obj.hSlider.t(i).HorizontalAlignment = 'left';
            end
        end
        
        % Update parameters to be used as model input
        function updateparameters(obj)
            obj.ParamValueList = zeros(obj.nParam,1);
            % set parameter values based on current slider values
            for i = 1:obj.nParam
                sval = obj.hSlider.s(i).Value;
                obj.ParamValueList(i,1) = sval;
                obj.Parameters.(obj.ParamNameList{i}) = sval;
            end
        end
        
        % Run model with current parameters
        function solvemodel(obj)
            % time span for simulation
            tspan = 0:600;
            % start of IP3 stimulation
            stimonset = 200;
            obj.StimOnset = stimonset;
            P = obj.Parameters;
            % run model with parameter stucture
            [t, c, e, m, u, fluxes] = CAmodel(P, 'Tspan', tspan,...
                'StimOnset', stimonset);
            % grab output and store into object
            obj.Time = t;
            obj.StateVar(:,1) = c;
            obj.StateVar(:,2) = e;
            obj.StateVar(:,3) = m;
            obj.StateVar(:,4) = u;
            obj.NonStateVar = fluxes;   
        end
        
        % Establish initial axis properties
        function setinitialaxis(obj)
            ax = obj.hAxis.h(4); plot(ax, 0, 0)
            obj.setproperties(ax)
            xlabel(ax, 'Time (s)')
            ylabel(ax, '[Ca^{2+}]_{Mt} (\muM)')
            ax = obj.hAxis.h(3); plot(ax, 0, 0)
            obj.setproperties(ax)
            xlabel(ax, 'Time (s)')
            ylabel(ax, '[Ca^{2+}]_{Cyt} (\muM)')
            ax = obj.hAxis.h(2); plot(ax, 0, 0)
            obj.setproperties(ax)
            xlabel(ax, 'Time (s)')
            ylabel(ax, '[Ca^{2+}]_{\mud} (\muM)')
            ax = obj.hAxis.h(1); plot(ax, 0, 0)
            obj.setproperties(ax)
            xlabel(ax, 'Time (s)')
            ylabel(ax, '[Ca^{2+}]_{ER} (\muM)')
        end
        
        % Update plots based on current model state
        function updateaxes(obj,~,~)
            % update things and solve model
            updateparameters(obj);
            updatetext(obj);
            solvemodel(obj);     
            % grab state variables
            t = obj.Time;
            c = obj.StateVar(:,1);
            e = obj.StateVar(:,2);
            m = obj.StateVar(:,3);
            u = obj.StateVar(:,4);
            % update X and Y data of axes
            ax = obj.hAxis.h(4);
            ax.Children.XData = t;
            ax.Children.YData = m;
            ax = obj.hAxis.h(3);
            ax.Children.XData = t;
            ax.Children.YData = c;
            ax = obj.hAxis.h(2);
            ax.Children.XData = t;
            ax.Children.YData = u;
            ax = obj.hAxis.h(1);
            ax.Children.XData = t;
            ax.Children.YData = e;
        end
        
        % Set button positions
        function setbuttonposition(obj)
            nCol = 3;
            nRow = ceil(obj.nButton/3);
            width = 1/nCol;
            height = 1/nRow;
            ibutton = 0;
            for irow = 1:nRow
                for icol = 1:nCol
                    ibutton = ibutton + 1;
                    if ibutton > obj.nButton
                        break
                    end
                    obj.hButton.b(ibutton).Position = ...
                        [(icol-1)*width (nRow-irow)*height...
                        width height];
                end
            end
            
        end 
    end
    
    methods (Static)
        % Set axis properties 
        function setproperties(ax)
            set(gcf, 'Color', 'w')
            ax.FontName = 'arial';
            ax.FontSize = 12;
            ax.LineWidth = 1.2;
            box(ax, 'off')
        end
    end
    
    %% CREATE BUTTONS
    methods
        % function for adding buttons
        function addButton(obj, name, callback)
            obj.nButton = obj.nButton + 1;
            obj.hButton.b(obj.nButton) = ...
                uicontrol('Parent', obj.hButtonPanel.Controls,...
                'Style', 'pushbutton',...
                'String', name,...
                'FontSize', 10,...
                'Units', 'normalized',...
                'Callback', callback);
        end
    end
    
    %% DEFINE BUTTON CALLBACKS
    methods
        %  resets axes to initial state
        function onReset(obj,~,~)
            arrayfun(@cla,findall(0,'type','axes'))
            obj.Parameters = obj.InitialParam;
            obj.ParamNameList = fieldnames(obj.Parameters);
            obj.ParamValueList = zeros(obj.nParam,1);
            obj.nParam = numel(obj.ParamNameList);
            for i = 1:obj.nParam
                obj.ParamValueList(i) = ...
                    obj.Parameters.(obj.ParamNameList{i});
            end
            setsliders(obj)
            setsliderposition(obj)
            setinitialaxis(obj)
            updateaxes(obj)
            
        end
        % export the current parameter set
        function onExport(obj,~,~)
            ParamSet = obj.Parameters;
            path = uigetdir(pwd,'Select path for output');
            name = char(inputdlg("Enter parameter set name"));
            save([path '\' name],'ParamSet')
        end
        % import a parameter set
        function onImport(obj, ~, ~)
            [pSet,path] = uigetfile;
            pset = load([path pSet]);
            fname = char(fieldnames(pset));
            obj.hFigure.Name = pSet;
            obj.Parameters = pset.(fname);
            for i = 1:obj.nParam
                obj.hSlider.s(i,1).Value = ...
                    obj.Parameters.(obj.ParamNameList{i});
            end
            updateaxes(obj)
        end
        
        % show non-state variables in subplot
        function onShowNonState(obj,~,~)
            NonState = obj.NonStateVar;
            t = obj.Time;
            figure
            fields = fieldnames(NonState);
            for i = 1:numel(fields)
                VarXaxis = t;
                VarYaxis = obj.NonStateVar.(fields{i});
                subplot(3,6,i)
                plot(VarXaxis,VarYaxis, 'Color', 'k');                
                set(gcf, 'Color', 'w')
                ax = gca;
                ax.FontName = 'arial';
                ax.FontSize = 10;
                ax.Title.String = fields{i};
                ax.LineWidth = 1;
                ax.XAxis.Label.String = 'Time (sec)';
                box off
            end
        end
        
        % show bifurcation diagram
        function onBifDiag(obj,~,~)
            P = obj.Parameters;
            nstim = 80;
            minIP3 = 0;
            maxIP3 = 3;
            ip3 = linspace(minIP3, maxIP3, nstim);
            cmin = zeros(nstim, 3); cmax = zeros(nstim, 3);
            f = waitbar(0,'Creating bifurcation diagram...');
            
            dist = P.distance*1e-9; % [m]
            r = 0.58*1e-6; % [m]
            SA = 2*0.2*4*pi*r^2; % [m^2]
            nMtObject = 200;
            volMd = SA*nMtObject*dist; % [m^3]
            volMd = volMd * 1e15; % [pL]            
            tspan = 0:0.1:600;
            
            for i = 1:nstim
                P.ip3 = ip3(i);
                [t, c, ~, ~, u, ~] = CAmodel(P, 'Tspan', tspan);
                c = (P.volCt*c + u*volMd)/(P.volCt + volMd);
                cmax(i,1) = max(c(t>(t(end) - 100)));
                cmin(i,1) = min(c(t>(t(end) - 100)));
                waitbar(i/nstim,f)
            end
            close(f)            
            figure
            hold on
            plot(ip3, cmax, 'LineWidth',1.2,'Color','k')
            plot(ip3, cmin, 'LineWidth',1.2,'Color','k')
            set(gcf, 'Color', 'w')
            ax = gca;
            ax.FontName = 'arial';
            ax.FontSize = 16;
            ax.LineWidth = 1.5;
            xlabel('[IP_3] (\muM)')
            ylabel('[Ca^2^+]_{Cyt^{eff}} (\muM)')
            box off
        end
        
        % callback for normalized change in calcium from basal levels
        function onNormalized(obj, ~, ~)
            % grab state variables 
            t = obj.Time;
            c = obj.StateVar(:,1);
            e = obj.StateVar(:,2);
            m = obj.StateVar(:,3);
            u = obj.StateVar(:,4);
            tend = t(end);
            timeWindow = 0.2*tend;
            tplot = t(t>(t(end) - timeWindow));
            ERPlot = e(t>(t(end) - timeWindow)); % cycle based on ER
            CytPlot = c(t>(t(end) - timeWindow));
            MtPlot = m(t>(t(end) - timeWindow));
            MdPlot = u(t>(t(end) - timeWindow));

            % get index right before stimulation
            tInd = find(t < obj.StimOnset, 1, 'last');
            % cytosol normalized change relative to basal
            cSS = c(tInd);
            deltaC = CytPlot - cSS;
            maxDeltaC = max(abs(deltaC));
            deltaCNORM = deltaC/maxDeltaC;
            % ER normalized change relative to basal
            eSS = e(tInd);
            deltaE = ERPlot - eSS;
            maxDeltaE = max(abs(deltaE));
            deltaENORM = 1 + deltaE/maxDeltaE;
            % Mt normalized change relative to basal
            mSS = m(tInd);
            deltaM = MtPlot - mSS;
            maxDeltaM = max(abs(deltaM));
            deltaMNORM = deltaM/maxDeltaM;
            % microdomain normalized change relative to basal
            uSS = u(tInd);
            deltaU = MdPlot - uSS;
            maxDeltaU = max(abs(deltaU));
            deltaUNORM = deltaU/maxDeltaU;
            % create figure
            figure('Color', 'w')
            % plot 
            ax = gca;
            h = plot(tplot,[deltaCNORM,...
                deltaENORM, deltaMNORM,...
                deltaUNORM], 'Parent',ax);
            % set line properties
            h(1).LineWidth = 1;
            h(2).LineWidth = 1;
            h(3).LineWidth = 1;
            h(4).LineWidth = 1;
            h(1).Color = [0, 0, 0];
            h(2).Color = [0.75, 0.0780, 0.1840];
            h(3).Color = [0.3010, 0.7450, 0.9330];
            h(4).Color = [0, 0.5, 0];
            % axis properties
            ax.FontName = 'arial';
            ax.FontSize = 18;
            ax.LineWidth = 1.2;
            box off
            xlabel('Time (s)')
            ylabel('Relative Change from Basal')
        end
        
        
    end
end

