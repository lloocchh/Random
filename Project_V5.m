%% Initialisation block 
close all;
clc;
clear;

%% Your information 
disp("Student name : Loch Sharkey");
disp("   Student ID: 220219776");
disp("Project Code For SEN771 T1 2023");
disp("X Factor 1: Random target locations and number of targets");
disp("X Factor 2: Multiplier for more precise pathing at cost of computer resources");
disp("Warning: Setting 500 target points or a multiplier above 8 creates lots of lag")
disp("-------------------------------");

%% Variables
Height = 90;
Width = 120;
Multipier = 2; % Minimum of 2

% Grid Options
defaultOption = '0';
gridEnabled = inputdlg('Enable grid? (1 / 0):', 'Grid Options', 1, {defaultOption});
gridEnabled = str2double(gridEnabled{1});
    
% Check if the input is either true or false
if ~(gridEnabled == 1 || gridEnabled == 0)
    error('Invalid input. Please enter 1 for true or 0 for false.');
end

% Obtacle Options
defaultOption = '10';
ObsWidth = inputdlg('Obstacle Width (0 - 10):', 'Grid Options', 1, {defaultOption});
ObsWidth = str2double(ObsWidth{1});
    
% Check if the input is within the valid range
if ObsWidth < 0 || ObsWidth > 10
    error('Invalid input. Obstacle width must be between 0 and 10 meters.');
end

ObsLen=2*ObsWidth; 

% Target Options
defaultOption = '1';
RandTargetEnabled = inputdlg('Random Target Locations? (1 / 0):', 'Grid Options', 1, {defaultOption});
RandTargetEnabled = str2double(RandTargetEnabled{1});
    
% Check if the input is empty (user pressed Enter) and set the default value
if isempty(RandTargetEnabled)
    RandTargetEnabled = 1;
else
    % Check if the input is either true or false
    if ~(RandTargetEnabled == 1 || RandTargetEnabled == 0)
        error('Invalid input. Please enter 1 for true or 0 for false.');
    end
end

if (RandTargetEnabled == 1)
    defaultOption = '10';
    NumTargets = inputdlg('Number of Target Points (1 - 500):', 'Grid Options', 1, {defaultOption});
    NumTargets = str2double(NumTargets{1});
end

% Program Variables
global TargetsFound;
TargetsFound = 0;
Offset = 0.75;

Connecting_Distance = 2; 

MAP = int8(zeros(Height * Multipier, Width * Multipier));

%% Scenario Field dimentions and visulisation
Xcorners=[0,Width, Width, 0,0];  % standard soccer field and they are im meters
Ycorners=[0,0,Height,Height, 0];
Fhandle = figure();
plot(Xcorners,Ycorners, 'LineWidth',3)
axis([-20 Width+20 -20 Height+20]) ;
xlabel("x(m)");
ylabel("y(m)");
title('َSoccer Field')
hold on

%% Obstacle locations 
O(:,1)=[30,70]'; 
O(:,2)=[30, 45]';
O(:,3)=[30, 20]';
O(:,4)=[45,60]';
O(:,5)=[45,30]';
O(:,6)=[60,80]' ;
O(:,7)=[60,65]';
O(:,8)=[60,25]';
O(:,9)=[60,10]';
O(:,10)=[75,69]'; 
O(:,11)=[75,30]';
O(:,12)=[90,70]';
O(:,13)=[90,45]';
O(:,14)=[90,20]';
O(:,15)=[109,45]';
plot(O(1,:),O(2,:),'d');

%% plot obstacles
for i=1:15
    if(i<6 || i>9)   % Thin rectangles
        rectangle('Position',[O(1,i)-ObsWidth/2 O(2,i)-ObsWidth, ObsWidth,ObsLen],'FaceColor',[1 0 0],'EdgeColor','r');
    else          % Fat rectangles
        rectangle('Position',[O(1,i)-ObsWidth O(2,i)-ObsWidth/2, ObsLen, ObsWidth],'FaceColor',[1 0 0],'EdgeColor','r');
    end
end

%% Add Obstructions To Map as Ones
for i=1:15
    if(i<6 || i>9)
        X1 = round((O(1,i) - ObsWidth/2 - Offset) * Multipier);
        X2 = round(X1 + ((ObsWidth + Offset*2) * Multipier));
        Y1 = round((O(2,i) - ObsWidth - Offset) * Multipier);
        Y2 = round(Y1 + ((ObsLen + Offset*2) * Multipier));
    else
        X1 = round((O(1,i) - ObsWidth - Offset) * Multipier);
        X2 = round(X1 + ((ObsLen + Offset*2) * Multipier));
        Y1 = round((O(2,i) - ObsWidth/2 - Offset) * Multipier);
        Y2 = round(Y1 + ((ObsWidth + Offset*2) * Multipier));
    end
    MAP(Y1:Y2, X1:X2) = 1;
end

%% Targets 
if (RandTargetEnabled == true)
    target_points = [];
    while length(target_points) < NumTargets
        % Generate random location for target point
        x = randi(Width - 2) + 1;
        y = randi(Height- 2) + 1;

        is_under_obstacle = false;
        % Check if the target point is under any obstacle or another target
        % point
        if ~(MAP(y * Multipier, x * Multipier) == 0)
            is_under_obstacle = true;
        end
    
        % If the target point is not under any obstacle, add it to the list
        if ~is_under_obstacle
            target_points = [target_points; x, y];
            Target_Plot = plot(x,y,'og','LineWidth',1);

            % Update the MAP to indicate the target point
            MAP(y * Multipier, x * Multipier) = 2;
        end
    end
    % Display the target points
    disp('Target Points:');
    disp(target_points);

else
    NumTargets = 4;
    % Plot the given target points
    Target(:,1)=[119.5,54.5]';
    Target(:,2)=[82.5, 15.5]';
    Target(:,3)=[82.5, 75.5]';
    Target(:,4)=[60.5, 54.5]';

    Target_Plot = plot(Target(1,:),Target(2,:),'og','LineWidth',1);

    % Add Goals To Map as Twos
    MAP(round(Target(2,1) * Multipier),round(Target(1,1) * Multipier)) = 2;
    MAP(round(Target(2,2) * Multipier),round(Target(1,2) * Multipier)) = 2;
    MAP(round(Target(2,3) * Multipier),round(Target(1,3) * Multipier)) = 2;
    MAP(round(Target(2,4) * Multipier),round(Target(1,4) * Multipier)) = 2;

end

%% Robot position
Robot=[0.5,45.5]';
Start_Plot = plot(Robot(1),Robot(2),'ok');

%% Show Grids (only for visualisation, remove later)
if (gridEnabled == true)   
     GridHLines=[zeros(Height+1,1) Width*ones(Height+1,1) [0:Height]' [0:Height]'];
     GridVlines=[[0:Width]' [0:Width]' zeros(Width+1,1) Height*ones(Width+1,1)];
     
     for i=1:Height
         line(GridHLines(i,[1 2]),GridHLines(i,[3 4]))
     end
     
     for i=1:Width
         line(GridVlines(i,[1 2]),GridVlines(i,[3 4]))
     end
end

%% MAP Creation
% Create Variables and Map Array
disp("      Number of Goals: " + string(NumTargets));
disp("      Grid Multiplier: " + string(Multipier));
disp("         Offset: " + string(Offset) + "m");
disp("-------------------------------");

%% Define the Robot Start Position
StartX = Robot(1) * Multipier;
StartY = Robot(2) * Multipier;

%% Run A* Algorithum To Find Goals
for i=1:NumTargets
    % Find path to closest goal
    OptimalPath = findOptimalPath(StartX, StartY, Height * Multipier, Width * Multipier, MAP, Connecting_Distance);
    
    % If no path to the next goal exists then give up and return home
    if CheckPath(OptimalPath) == 1
        break
    end
    
    % Set the goal as the new start
    StartX = OptimalPath(1,2);
    StartY = OptimalPath(1,1);
    
    % Remove the goal from the list
    MAP(OptimalPath(1,1),OptimalPath(1,2)) = 3;

    % Check if this is the first iteration
    if exist('TotalPath') == 1
        TotalPath = vertcat(OptimalPath, TotalPath);
    else
        TotalPath = OptimalPath;    
    end
end

%% Run A* Algorithum To Complete Loop
MAP(45.5 * Multipier, 0.5 * Multipier) = 2;
OptimalPath = findOptimalPath(StartX, StartY, Height * Multipier, Width * Multipier, MAP, Connecting_Distance);
TotalPath = vertcat(OptimalPath, TotalPath);

TotalPath_Plot = plot(TotalPath(:,2)./Multipier,TotalPath(:,1)./Multipier,'.-r','color','k');

%% Calculate The Total Path Length
d = diff(TotalPath); % EDIT
TotalPath_Length = round(sum(sqrt(sum(d.*d,2))/Multipier),2);

disp("       Path Length: " + string(TotalPath_Length) + "m");

%% Time and time step values
TimeIndex = 1;
dT = 1;  % reasonable choice of dT, for the environment that things dont move fast
TotalPath_Points = size(TotalPath, 1);

%% Simulation Runtime Visualization
for K = TotalPath_Points:-dT:1
    TimeHandler = text(100, 100, ["Time (s) : "+string(TotalPath_Points - K)]);
    TimeIndex = K;
    text
    
    % Check if this is the first iteration. If so crate the robot
    if K == TotalPath_Points
        RobotLoc(:,TimeIndex) = Robot;
        Roh = plot(RobotLoc(1,TimeIndex),RobotLoc(2,TimeIndex),'om');
        legend([Start_Plot Roh Target_Plot TotalPath_Plot],'Start', 'Robot', 'Goal', 'Path','Location','northwest')

    % Otherwise, delete the robot and move it to the next coodinate
    else
        RobotLoc(:,TimeIndex) = TotalPath(TimeIndex,:)';
        set(Roh,'XData',RobotLoc(2,TimeIndex)./Multipier,'YData',RobotLoc(1,TimeIndex)./Multipier);
    end
    
    drawnow
    pause(0.01);
    delete(TimeHandler); 
end

disp("      Total Time Taken: " + string(TotalPath_Points) + "s");
disp("-------------------------------");

%% Functions
% Check if the algorithum failed to find a path to the next goal
function output = CheckPath(input)
    global TargetsFound;
    if input(1) == inf
        h=msgbox('Sorry, No path exists to the Target!','warn');
        output = 1;
    else
        TargetsFound = TargetsFound + 1;
        output = 0;
    end
end
function OptimalPath = findOptimalPath(startX, startY, height, width, map, connectingDistance)
    % A* Pathfinding Algorithm
    % Inputs:
    %   startX: Starting X position
    %   startY: Starting Y position
    %   map: A matrix representing the map.
    %        0 represents a passable cell
    %        1 represents an obstacle
    %        2 represents a target point
    %   connectingDistance: Maximum distance to consider neighboring cells
    
    % Output:
    %   OptimalPath: An array of (x, y) tuples representing the optimal path
    %   from the starting position to the goal, or an empty array if no path exists.

    % Variable Setup
    gScore = zeros(height, width);
    fScore = single(inf(height, width));
    hn = single(zeros(height, width));
    openMat = int8(zeros(height, width));
    closedMat = int8(zeros(height, width));
    closedMat(map == 1) = 1;
    parentX = int16(zeros(height, width));
    parentY = int16(zeros(height, width));

    % Set up matrices representing neighbors to be investigated
    neighborCheck = ones(2 * connectingDistance + 1);
    dummy = 2 * connectingDistance + 2;
    mid = connectingDistance + 1;
    for i = 1:connectingDistance - 1
        neighborCheck(i, i) = 0;
        neighborCheck(dummy - i, i) = 0;
        neighborCheck(i, dummy - i) = 0;
        neighborCheck(dummy - i, dummy - i) = 0;
        neighborCheck(mid, i) = 0;
        neighborCheck(mid, dummy - i) = 0;
        neighborCheck(i, mid) = 0;
        neighborCheck(dummy - i, mid) = 0;
    end
    neighborCheck(mid, mid) = 0;

    [row, col] = find(neighborCheck == 1);
    neighbors = [row col] - (connectingDistance + 1);
    numNeighbors = size(col, 1);

    % Create heuristic matrix based on distance to nearest goal node
    [col, row] = find(map == 2);
    registeredGoals = [row col];
    numGoals = size(registeredGoals, 1);

    for k = 1:height
        for j = 1:width
            if map(k, j) == 0
                mat = registeredGoals - repmat([j k], numGoals, 1);
                distance = min(sqrt(sum(abs(mat).^2, 2)));
                hn(k, j) = distance;
            end
        end
    end

    % Initialize start node with fScore and open first node
    fScore(startY, startX) = hn(startY, startX);
    openMat(startY, startX) = 1;

    while true
        minOpenFScore = min(min(fScore));
        if minOpenFScore == inf
            % Failure!
            OptimalPath = [inf];
            reconstructPath = false;
            break;
        end
        [currentY, currentX] = find(fScore == minOpenFScore);
        currentY = currentY(1);
        currentX = currentX(1);

        if map(currentY, currentX) == 2
            % Goal found!
            reconstructPath = true;
            break;
        end

        % Remove node from open list to closed list
        openMat(currentY, currentX) = 0;
        fScore(currentY, currentX) = inf;
        closedMat(currentY, currentX) = 1;
        for p = 1:numNeighbors
            i = neighbors(p, 1); % Y
            j = neighbors(p, 2); % X
            if currentY + i < 1 || currentY + i > height || currentX + j < 1 || currentX + j > width
                continue
            end
            flag = 1;
            if closedMat(currentY + i, currentX + j) == 0 % Neighbor is open
                if abs(i) > 1 || abs(j) > 1
                    % Need to check that the path does not pass an object
                    jumpCells = 2 * max(abs(i), abs(j)) - 1;
                    for k = 1:jumpCells
                        yPos = round(k * i / jumpCells);
                        xPos = round(k * j / jumpCells);

                        if map(currentY + yPos, currentX + xPos) == 1
                            flag = 0;
                        end
                    end
                end

                if flag == 1
                    tentativeGScore = gScore(currentY, currentX) + sqrt(i^2 + j^2);
                    if openMat(currentY + i, currentX + j) == 0
                        openMat(currentY + i, currentX + j) = 1;
                    elseif tentativeGScore >= gScore(currentY + i, currentX + j)
                        continue
                    end
                    parentX(currentY + i, currentX + j) = currentX;
                    parentY(currentY + i, currentX + j) = currentY;
                    gScore(currentY + i, currentX + j) = tentativeGScore;
                    fScore(currentY + i, currentX + j) = tentativeGScore + hn(currentY + i, currentX + j);
                end
            end
        end
    end

    k = 2;
    if reconstructPath
        optimalPath(1, :) = [currentY currentX];
        while reconstructPath
            currentXDummy = parentX(currentY, currentX);
            currentY = parentY(currentY, currentX);
            currentX = currentXDummy;
            optimalPath(k, :) = [currentY currentX];
            k = k + 1;
            if currentX == startX && currentY == startY
                break
            end
        end
    end
    
    OptimalPath = optimalPath;
end
