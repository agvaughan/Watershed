function h_out = namedFigure(varargin)
%{ 
    Searches for as figure with the given name; 
    If it doesn't find one, it creates a new figure.
%}

% Treat inputs as inputs to sprintf.
thisName = sprintf(varargin{:});

%All figures
figureList = allchild(0);

%Figure handle of all figures matchin thisName
figuresWithThisName = sort(figureList(strcmp(get(figureList,'Name'),thisName)));

%% Load or switch to the appropriate figure.
switch length(figuresWithThisName)
    case 0,
        % None found; make a new figure
        h = figure('Name',thisName);
        set(gca,'visible', 'off');
        
        if strfind(thisName,'Behavioral Response')
            %set(this_gcf, 'PaperSize', [10 15]);
        else
            %set(h, 'color', [1 1 1]);
            set(h,'color','w');
        end
        
        %axis square;

    case 1
        % Switch to that figure.
        figure(figuresWithThisName);
        h = figuresWithThisName;
    otherwise
        % Over-write the last figure, I guess?
        figure(figuresWithThisName(end));
        h = figuresWithThisName(end);
end

if nargout,
    h_out = h;
end


