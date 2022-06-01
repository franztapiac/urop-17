%% Processing of Force and XYZ DATA from Vicon Nexus

% Aim:
% To compare the efficiency in motion between an average healthy subject
% and an average patient with a prosthetic

% Objectives:
% 1. Compare average distance covered for movement between subjects
    % Hypothesis: For subjects with prosthetics, greater distance was covered
    %   when performing trial task.
    % Plot 1.1: Bar chart comparing the average, w/ error bars
    % Plot 1.2: Average distance covered VS Time elapsed, but instead of plotting a
    %   line, plot a point that moves as time increases. (may help to see
    %   when in time, or when during the movement, the distance largely changes
% 2. Compare average sum of angles swept for movement between subjects
    % Hypothesis: Healthy subjects will have required minimal joint readjustments
    % Plot 2.1: Bar chart
    % Plot 2.2: Angle VS Time, as in Plot 1.2
% 3. Compare average energy expenditure between subjects
    % Plot 3.1: Bar chart
    % Plot 3.2: Energy VS Time. Like in Plot 1.2
% 4. Compare average changes in posture between subjects
    % Plot 4.1: X-Y graph showing oscillations in posture (COM) during
    % trial.

% Processing Rate details
frameRate = 200; % frame/s
sampleRate = 1000; % Hz, 1/s

% Details for Low-Pass Filter
filterOrder = 3;        % 60 Db/dec attenuation after cut-off frequency 
cutOffFrequency = 3;    % rad/s
[numerator,denominator] = butter(filterOrder, cutOffFrequency/(frameRate/2), 'low');    % Results in Transfer Function coefficients of filter

% hfvt = fvtool(numerator,denominator,'FrequencyScale','log'); % Plots the Mag Response of the filter

% Trial data files
taskString1 = {'CPRT','CPRT 1','CPRT 2'};
taskString2 = {'BandB','BandB 1','BandB 2'};
taskString3 = {'light_power','heavy_power','empty_tin','full_jar','door_handle','screw', 'food_cutting'};
taskString4 = {'light_tripod','heavy_tripod', 'button_board','food_cutting'};
taskString5 = {'light_sphere','heavy_sphere','jar_lid','carton_pour'};
taskString6 = {'light_lateral','heavy_lateral','jug_pour','key','zip','tray'};
taskString7 = {'light_tip','heavy_tip','coins','zip','button_board'};
taskString8 = {'light_extension','heavy_extension','page_turning','tray'};

taskString = {taskString1, taskString2, taskString3, taskString4, taskString5, taskString6, taskString7, taskString8};

% Global variables used throughout the code
subChar_cell = {'H','S'};
subCount_cell = {[4,6:11],[7,13]};
subSesh_cell = {1,[1,2,3]};
trials_cell = {1,3};

% For loops for subject and trial data files
groupCounter = 1;
for subType = 1:2
    subChar = subChar_cell{subType};
    subCount = subCount_cell{subType};
    subSesh = subSesh_cell{subType};
    trials = trials_cell{subType};

    for subID = subCount
        for session = subSesh
            for taskID = 1:trials
                % Accessing spreadSheet and extracting data
                if subID >= 10
                    if subType == 1
                        spreadSheet = [cd,'\data\' subChar num2str(subID) '\Session ' num2str(session) '\' taskString{groupCounter}{taskID} '.csv'];
                    else
                        spreadSheet = [cd,'\data\' subChar num2str(subID) '\Session ' num2str(session) '\' taskString{groupCounter}{taskID} '.csv'];
                    end
                else
                    if subType == 1
                        spreadSheet = [cd,'\data\' subChar '0' num2str(subID) '\Session ' num2str(session) '\' taskString{groupCounter}{taskID} '.csv'];
                    else
                        spreadSheet = [cd,'\data\' subChar '0' num2str(subID) '\Session ' num2str(session) '\' taskString{groupCounter}{taskID} '.csv'];
                    end
                end
                
                % Load all the data onto MatLab
                [numData, txtData, ~] = xlsread(spreadSheet);
                
                % XYZ DATA in millimeters, mm
                tmpTime = numData(:, 1);                                    % temporary data set
                xyzStart = find(isnan(tmpTime),1,'last');                   % 1 rows above where XYZ DATA starts ON the Matrix copy of Excel spreadsheet
                frameStart = tmpTime(xyzStart+1,1);                         % number of starting frame, may not be 1
                frameEnd = tmpTime(end);                                	% number of ending frame
                clear tmpTime;
                xyzData = numData(xyzStart+1:end,3:53);                     

                % XYZ DATA labels
                tmpLabels = txtData(xyzStart-1,3:end-2);
                for iter = 1:3:size(tmpLabels,2)
                    [~,tmp] = strread(tmpLabels{iter}, '%s %s', 'delimiter',':');   % e.g. read 'Patient30' & 'LACR', separated by ':', but discard 'Patient30'
                    markerLabels{floor(iter./3)+1} = tmp{1};                         % isolate the marker names by rounded divisio
                end

                % XYZ DATA Filtering
                xyzData_f = filtfilt(numerator,denominator,xyzData);    % Apply butterworth filter to XYZ data to remove noise
                clear xyzData txtData;

                % Defining the Observed Side
                defineSide = find(ismember(markerLabels,'LPALML'));
                if isempty(defineSide)
                   observedSide = 'R';
                else
                   observedSide = 'L';
                end

                % Find indeces of the corresponding labels
                % PALM
                palmInd(1) = find(ismember(markerLabels,[observedSide 'PALMM'])); % observedSide allows for a general code for right or left side
                palmInd(2) = find(ismember(markerLabels,[observedSide 'PALMC']));
                palmInd(3) = find(ismember(markerLabels,[observedSide 'PALML']));

                % FOREARM
                forearmInd(1) = find(ismember(markerLabels,[observedSide 'FARML']));
                forearmInd(2) = find(ismember(markerLabels,[observedSide 'FARMM']));
                forearmInd(3) = find(ismember(markerLabels,[observedSide 'FARMC']));

                % UPPERARM
                upperarmInd(1) = find(ismember(markerLabels,[observedSide 'ARM1']));
                upperarmInd(2) = find(ismember(markerLabels,[observedSide 'ARM2']));
                upperarmInd(3) = find(ismember(markerLabels,[observedSide 'ARM3']));

                % TRUNK
                trunkInd(1) = find(ismember(markerLabels,'LACR'));
                trunkInd(2) = find(ismember(markerLabels,'RACR'));
                trunkInd(3) = find(ismember(markerLabels,'STER'));
                trunkInd(4) = find(ismember(markerLabels,'C7'));

                % ELBOW
                elbowInd(1) = find(ismember(markerLabels,[observedSide 'EPICEM']));
                elbowInd(2) = find(ismember(markerLabels,[observedSide 'EPICEL']));

                % WRIST
                wristInd(1) = find(ismember(markerLabels,[observedSide 'EPICWM']));
                wristInd(2) = find(ismember(markerLabels,[observedSide 'EPICWL']));

                % Go through all the frames and find the centroids (centre of location/mass)
                for iter = 1:size(xyzData_f,1)     % iterate through all (remaining - if cut) frames

                    % these calculations rely on the fact that there are no gaps, or any NaNs
                    palmCentroid(iter,:) = [...
                        (xyzData_f(iter,(palmInd(1)-1)*3+1) + xyzData_f(iter,(palmInd(2)-1)*3+1) + xyzData_f(iter,(palmInd(3)-1)*3+1))/3, ...
                        (xyzData_f(iter,(palmInd(1)-1)*3+2) + xyzData_f(iter,(palmInd(2)-1)*3+2) + xyzData_f(iter,(palmInd(3)-1)*3+2))/3, ...
                        (xyzData_f(iter,(palmInd(1)-1)*3+3) + xyzData_f(iter,(palmInd(2)-1)*3+3) + xyzData_f(iter,(palmInd(3)-1)*3+3))/3  ...
                                           ];      % xyzData_f contains XYZ data of each marker in each next 3 columns
                                                   % e.g. (palmInd(1)-1)*3+2 reps. the y-coordinate data of RPALMM (marker 15) on cell (1,44)

                    forearmCentroid(iter,:) = [...
                        (xyzData_f(iter,(forearmInd(1)-1)*3+1) + xyzData_f(iter,(forearmInd(2)-1)*3+1) + xyzData_f(iter,(forearmInd(3)-1)*3+1))/3, ...
                        (xyzData_f(iter,(forearmInd(1)-1)*3+2) + xyzData_f(iter,(forearmInd(2)-1)*3+2) + xyzData_f(iter,(forearmInd(3)-1)*3+2))/3, ...
                        (xyzData_f(iter,(forearmInd(1)-1)*3+3) + xyzData_f(iter,(forearmInd(2)-1)*3+3) + xyzData_f(iter,(forearmInd(3)-1)*3+3))/3  ...
                                              ];

                    upperarmCentroid(iter,:) = [...
                        (xyzData_f(iter,(upperarmInd(1)-1)*3+1) + xyzData_f(iter,(upperarmInd(2)-1)*3+1) + xyzData_f(iter,(upperarmInd(3)-1)*3+1))/3, ...
                        (xyzData_f(iter,(upperarmInd(1)-1)*3+2) + xyzData_f(iter,(upperarmInd(2)-1)*3+2) + xyzData_f(iter,(upperarmInd(3)-1)*3+2))/3, ...
                        (xyzData_f(iter,(upperarmInd(1)-1)*3+3) + xyzData_f(iter,(upperarmInd(2)-1)*3+3) + xyzData_f(iter,(upperarmInd(3)-1)*3+3))/3  ...
                                               ];

                    trunkCentroid(iter,:) = [...
                        (xyzData_f(iter,(trunkInd(1)-1)*3+1) + xyzData_f(iter,(trunkInd(2)-1)*3+1) + ...
                            xyzData_f(iter,(trunkInd(3)-1)*3+1) + xyzData_f(iter,(trunkInd(4)-1)*3+1))/4, ...
                        (xyzData_f(iter,(trunkInd(1)-1)*3+2) + xyzData_f(iter,(trunkInd(2)-1)*3+2) + ...
                            xyzData_f(iter,(trunkInd(3)-1)*3+2) + xyzData_f(iter,(trunkInd(4)-1)*3+2))/4, ...
                        (xyzData_f(iter,(trunkInd(1)-1)*3+3) + xyzData_f(iter,(trunkInd(2)-1)*3+3) + ...
                            xyzData_f(iter,(trunkInd(3)-1)*3+3) + xyzData_f(iter,(trunkInd(4)-1)*3+3))/4 ...
                                            ];

                    wristCentroid(iter,:) = [...
                        (xyzData_f(iter,(wristInd(1)-1)*3+1) + xyzData_f(iter,(wristInd(2)-1)*3+1))/2, ...
                        (xyzData_f(iter,(wristInd(1)-1)*3+2) + xyzData_f(iter,(wristInd(2)-1)*3+2))/2, ...
                        (xyzData_f(iter,(wristInd(1)-1)*3+3) + xyzData_f(iter,(wristInd(2)-1)*3+3))/2  ...
                                            ];
                                        
                    if observedSide == 'R'
                        flag = 4;
                    else
                        flag = 3;
                    end
                    shoulderAngle(iter,:)=angle_among_three_points(...
                        trunkCentroid(iter,:),...
                        xyzData_f(iter,(trunkInd(flag)-1)*3+1:(trunkInd(flag)-1)*3+3),...
                        xyzData_f(iter,(elbowInd(1)-1)*3+1:(elbowInd(1)-1)*3+3));
                    elbowAngle(iter,:)=angle_among_three_points(...
                        xyzData_f(iter,(trunkInd(flag)-1)*3+1:(trunkInd(flag)-1)*3+3),...
                        xyzData_f(iter,(elbowInd(1)-1)*3+1:(elbowInd(1)-1)*3+3),...
                        wristCentroid(iter,:));
                                        
                end
                clear xyzData_f;

                % Removing the starting bias from the centroids,
                % for each trial in each session, of each subject, of each type
                palmCentroidXYZ = palmCentroid-(repmat(mean(palmCentroid(1:50,:))',1,size(palmCentroid,1)))';
                %forearmCentroidXYZ{subType,subID,session,taskID} = forearmCentroid-(repmat(mean(forearmCentroid(1:50,:))',1,size(forearmCentroid,1)))';
                %upperarmCentroidXYZ{subType,subID,session,taskID} = upperarmCentroid-(repmat(mean(upperarmCentroid(1:50,:))',1,size(upperarmCentroid,1)))';
                %trunkCentroidXYZ{subType,subID,session,taskID} = trunkCentroid-(repmat(mean(trunkCentroid(1:50,:))',1,size(trunkCentroid,1)))';
                %wristCentroidXYZ{subType,subID,session,taskID} = wristCentroid-(repmat(mean(wristCentroid(1:50,:))',1,size(wristCentroid,1)))';
                %clear palmCentroid forearmCentroid upperarmCentroid trunkCentroid wristCentroid;

                % Radial displacement of centroids from their origin
                % for each trial in each session, of each subject, of each type
                palmCentroidDispPerSub = sqrt(sum(palmCentroidXYZ.^2,2));
                %forearmCentroidDispPerSub{subType,subID,session,taskID} = sqrt(sum(forearmCentroidXYZ{subType,subID,session,taskID}.^2,2));
                %upperarmCentroidDispPerSub{subType,subID,session,taskID} = sqrt(sum(upperarmCentroidXYZ{subType,subID,session,taskID}.^2,2));
                %trunkCentroidDispPerSub{subType,subID,session,taskID} = sqrt(sum(trunkCentroidXYZ{subType,subID,session,taskID}.^2,2));
                %wristCentroidDispPerSub{subType,subID,session,taskID} = sqrt(sum(wristCentroidXYZ{subType,subID,session,taskID}.^2,2));

                % Velocity of Centroids over each frame
                %palmCentroidAbsVlctPerSub = abs(diff([palmCentroidDispPerSub(1,:); palmCentroidDispPerSub])/(1/200));
%!!             All markers need to be complete, i.e. w/out gaps

                %Norm Start & End Frame calculation - make separate function
                %{

                % Calculation of start and end frames of normalization
                for frame = 1:size(xyzTime{subType,subID,session,taskID}(:,1))
                    if (palmCentroidDispPerSub{subType,subID,session,taskID}(frame,1) >= 10) && (palmCentroidAbsVlctPerSub{subType,subID,session,taskID}(frame,1) >= 39.6)
                        NormStartFrame{subType,subID,session,taskID} = frame;
                        break;
                    end
                end

                    % First, find the index of the Max point on Displacement graph
                [~,NormEndFrame{subType,subID,session,taskID}] = max(palmCentroidDispPerSub{subType,subID,session,taskID});   % the index is returned to: [~,Index] 
                    % Then, after that point of Max Displacement, find the point of
                    % Max Abs Velocity -> I want this one because this is the after the 4th clothespin is placed
                [~,AbsVlctMaxIndx] = max(palmCentroidAbsVlctPerSub{subType,subID,session,taskID}(NormEndFrame{subType,subID,session,taskID}:size(palmCentroidAbsVlctPerSub{subType,subID,session,taskID},1),1));
                    % Then because I made a smaller matrix ^to find the AbsVlctMaxIndx,
                    % I calculate the actual index again by:
                NormEndFrame{subType,subID,session,taskID} = NormEndFrame{subType,subID,session,taskID} + AbsVlctMaxIndx - 1;
                    % Now that I have an index after the 4th max peak and before returning
                    % to stationary position, from this point forward, finds the FIRST 
                    % occuring trough or if there is none, take the last frame number
                for frame = NormEndFrame{subType,subID,session,taskID}+1:size(palmCentroidDispPerSub{subType,subID,session,taskID},1)
                    if (palmCentroidDispPerSub{subType,subID,session,taskID}(frame,1) > palmCentroidDispPerSub{subType,subID,session,taskID}(frame-1,1)) && (palmCentroidDispPerSub{subType,subID,session,taskID}(frame-2,1) > palmCentroidDispPerSub{subType,subID,session,taskID}(frame-1,1))
                        NormEndFrame{subType,subID,session,taskID} = frame-1;
                        break;
                    elseif frame == size(palmCentroidDispPerSub{subType,subID,session,taskID},1)   % if frame gets to the end without finding a trough
                        NormEndFrame{subType,subID,session,taskID} = frame;
                        break;
                    end
                end
                %}

                %                                                                 <-First Dimension->
                %                       <-------------------------------------------Second dimension------------------------------------>
                %                       <----------------------- H ------------------------><-------------------- S -------------------->
                % Fourth Dimension: CPRT trial
                NormStartFrame        = {0, 0, 0, 542, 0, 679, 603, 881, 388, 570, 610, 0, 0; 0, 0, 0, 0, 0, 0, 159, 0, 0, 0, 0, 0,  526};  % Session 1, ^ Third
                NormStartFrame(:,:,2) = {0, 0, 0,   0, 0,   0,   0,   0,   0,   0,   0, 0, 0; 0, 0, 0, 0, 0, 0, 200, 0, 0, 0, 0, 0, 1489};  % Session 2, ^ dimension
                NormStartFrame(:,:,3) = {0, 0, 0,   0, 0,   0,   0,   0,   0,   0,   0, 0, 0; 0, 0, 0, 0, 0, 0, 733, 0, 0, 0, 0, 0,  591};  % Session 3, ^
                % ^ CPRT
                
                NormStartFrame(:,:,1,2) = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0,  72, 0, 0, 0, 0, 0,  649};
                NormStartFrame(:,:,2,2) = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 130, 0, 0, 0, 0, 0, 1042};
                NormStartFrame(:,:,3,2) = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 635, 0, 0, 0, 0, 0,  485};
                % ^ CPRT 1
                
                NormStartFrame(:,:,1,3) = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 547, 0, 0, 0, 0, 0, 516};
                NormStartFrame(:,:,2,3) = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0, 913};
                NormStartFrame(:,:,3,3) = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0,  80, 0, 0, 0, 0, 0, 494};
                % ^ CPRT 2
                
                NormEndFrame        = {0, 0, 0, 1427, 0, 1528, 1720, 1931, 1470, 1651, 1603, 0, 0; 0, 0, 0, 0, 0, 0,  5824, 0, 0, 0, 0, 0,  6013};
                NormEndFrame(:,:,2) = {0, 0, 0,    0, 0,    0,    0,    0,    0,    0,    0, 0, 0; 0, 0, 0, 0, 0, 0, 22607, 0, 0, 0, 0, 0, 12750};
                NormEndFrame(:,:,3) = {0, 0, 0,    0, 0,    0,    0,    0,    0,    0,    0, 0, 0; 0, 0, 0, 0, 0, 0, 10214, 0, 0, 0, 0, 0,  7923};
                
                NormEndFrame(:,:,1,2) = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0,  3443, 0, 0, 0, 0, 0,  7075};
                NormEndFrame(:,:,2,2) = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 13769, 0, 0, 0, 0, 0, 16841};
                NormEndFrame(:,:,3,2) = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0,  7798, 0, 0, 0, 0, 0,  6664};
                
                NormEndFrame(:,:,1,3) = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 3761, 0, 0, 0, 0, 0,  5386};
                NormEndFrame(:,:,2,3) = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 5931, 0, 0, 0, 0, 0, 15498};
                NormEndFrame(:,:,3,3) = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 7535, 0, 0, 0, 0, 0,  6643};
                
                % Cutting XYZ and Disp data
                palmCentroidXYZ_cut = palmCentroidXYZ(NormStartFrame{subType,subID,session,taskID}-frameStart+1:NormEndFrame{subType,subID,session,taskID}-frameStart+1,:);
                palmCentroidDispPerSub_cut = palmCentroidDispPerSub(NormStartFrame{subType,subID,session,taskID}-frameStart+1:NormEndFrame{subType,subID,session,taskID}-frameStart+1,:);
                %clear palmCentroidXYZ palmCentroidDispPerSub;
                
                palmCentroidXYZ_cut = palmCentroidXYZ_cut - (repmat(palmCentroidXYZ_cut(1,:)',1,size(palmCentroidXYZ_cut,1)))';
                palmCentroidDispPerSub_cut = palmCentroidDispPerSub_cut - palmCentroidDispPerSub_cut(1,:);
                
                palmCentroidXYZ_norm{subType,subID,session,taskID}= filtfilt(numerator,denominator,interp1(1:size(palmCentroidXYZ_cut,1),...
                                                                           palmCentroidXYZ_cut,...
                                                                           linspace(1,size(palmCentroidXYZ_cut,1),10000)))';              
                
                palmCentroidDispPerSub_norm{subType,subID,session,taskID}= filtfilt(numerator,denominator,interp1(1:size(palmCentroidDispPerSub_cut,1),...
                                                                           palmCentroidDispPerSub_cut,...
                                                                           linspace(1,size(palmCentroidDispPerSub_cut,1),10000)))';                

                % Calculating Angle of Shoulder and Elbow spanned during trial
                shoulderAngle_cut = shoulderAngle(NormStartFrame{subType,subID,session,taskID}-frameStart+1:NormEndFrame{subType,subID,session,taskID}-frameStart+1);
                elbowAngle_cut = elbowAngle(NormStartFrame{subType,subID,session,taskID}-frameStart+1:NormEndFrame{subType,subID,session,taskID}-frameStart+1);

                shoulderAngleArea{subType,subID,session,taskID} = sum(shoulderAngle_cut);
                elbowAngleArea{subType,subID,session,taskID} = sum(elbowAngle_cut);
                                                                       
                % Distance travelled by markers over each frame
                palmDist{subType,subID,session,taskID} = sum(euclidianDistanceCalculation(palmCentroidXYZ_cut));
                %forearmDist{subType,subID,session,taskID} = sum(euclidianDistanceCalculation(forearmCentroidDispPerSub{subType,subID,session,taskID}));
                %upperarmDist{subType,subID,session,taskID} = sum(euclidianDistanceCalculation(upperarmCentroidDispPerSub{subType,subID,session,taskID}));
                %trunkDist{subType,subID,session,taskID} = sum(euclidianDistanceCalculation(trunkCentroidDispPerSub{subType,subID,session,taskID}));
                clear palmCentroidXYZ_cut palmCentroidDispPerSub_cut;
                
                % Centre of Pressure (CoP) DATA in millimeters, mm
                if groupCounter == 1 || groupCounter == 2   % Standing trial: CPRT or BandB
                    
                    % Find out which force plate was used
                    forceData1_z = numData(5:54,5);
                    forceData2_z = numData(5:54,20);

                    if abs(mean(forceData1_z)) > abs(mean(forceData2_z))
                        copData_x = numData(5:(4+((frameEnd-frameStart+1)*sampleRate)/frameRate),9);
                        copData_y = numData(5:(4+((frameEnd-frameStart+1)*sampleRate)/frameRate),10);
                    else
                        copData_x = numData(5:(4+((frameEnd-frameStart+1)*sampleRate)/frameRate),24);
                        copData_y = numData(5:(4+((frameEnd-frameStart+1)*sampleRate)/frameRate),25);
                    end     % coordinates of force application point, from centre of force plate coordinate system
                    clear numData forceData1_z forceData2_z;

                    % CoP DATA Filtering
                    copData_x_f = filtfilt(numerator,denominator,copData_x);
                    clear copData_x;
                    copData_y_f = filtfilt(numerator,denominator,copData_y);
                    clear copData_y;
                    
                    if observedSide == 'R'
                        copData_x_f = -copData_x_f;     % When considering left-hand-performed trials
                        copData_y_f = -copData_y_f;
%!!                     ^ Comment out accordingly.
                    else
                        % copData_x_f = -copData_x_f;   % When considering right-hand-performed trials
                        % copData_y_f = -copData_y_f;
                    end
                    
                    % CoP DATA Normalization        
                    % Normalization to convert frame/time to percentage completion of task
                    copData_x_f_cut = copData_x_f(5*(NormStartFrame{subType,subID,session,taskID}-frameStart+1)-4:5*(NormEndFrame{subType,subID,session,taskID}-frameStart+1),:);
                    copData_y_f_cut = copData_y_f(5*(NormStartFrame{subType,subID,session,taskID}-frameStart+1)-4:5*(NormEndFrame{subType,subID,session,taskID}-frameStart+1),:);
                    
%??                 Would removing shift bias (see Static Standing CoP average) have a positive effect?
                    copData_x_filt{subType,subID,session,taskID} = copData_x_f_cut - copData_x_f_cut(1,:);
                    copData_y_filt{subType,subID,session,taskID} = copData_y_f_cut - copData_y_f_cut(1,:);
                    copData_disp_filt{subType,subID,session,taskID} = sqrt(copData_x_filt{subType,subID,session,taskID}.^2 + copData_y_filt{subType,subID,session,taskID}.^2);

                    copData_x_norm{subType,subID,session,taskID} = filtfilt(numerator,denominator,interp1(1:size(copData_x_filt{subType,subID,session,taskID},1),...
                                                                            copData_x_filt{subType,subID,session,taskID},...
                                                                            linspace(1,size(copData_x_filt{subType,subID,session,taskID},1),10000)))'; % 4000 = 5*800 frames
                    copData_y_norm{subType,subID,session,taskID} = filtfilt(numerator,denominator,interp1(1:size(copData_y_filt{subType,subID,session,taskID},1),...
                                                                            copData_y_filt{subType,subID,session,taskID},...
                                                                            linspace(1,size(copData_y_filt{subType,subID,session,taskID},1),10000)))';
                    copData_disp_norm{subType,subID,session,taskID} = filtfilt(numerator,denominator,interp1(1:size(copData_disp_filt{subType,subID,session,taskID},1),...
                                                                               copData_disp_filt{subType,subID,session,taskID},...
                                                                               linspace(1,size(copData_disp_filt{subType,subID,session,taskID},1),10000)))';
                    clear copData_x_f copData_y_f copData_x_f_cut copData_y_f_cut;
                end % if trial = CPRT or BandB
            end %taskID
        end %session
    end %subID
end % subType

%% Bar graphs of PalmCentroid Distance covered

for subType = 1:2
    if subType == 1
        palmDist_H = [palmDist{subType,:,1,1}]; % all values for able-bodied subject, Session 1
        palmDist_H(palmDist_H == 0) = [];
        palmDist_avg_H = mean(palmDist_H);
        palmDist_std_H = std(palmDist_H);
    else
        for subID = subCount_cell{subType}
            for session = subSesh_cell{subType}
                palmDist_avg_S{subID,session} = mean([palmDist{subType,subID,session,:}]);
                palmDist_std_S{subID,session} = std([palmDist{subType,subID,session,:}]);
            end %session
        end %subID
    end % if
end % subType

for subID = subCount_cell{2}
    palmDist_avg = [palmDist_avg_H palmDist_avg_S{subID,1};  palmDist_avg_H palmDist_avg_S{subID,2};  palmDist_avg_H palmDist_avg_S{subID,3}];
    palmDist_std = [palmDist_std_H palmDist_std_S{subID,1};  palmDist_std_H palmDist_std_S{subID,2};  palmDist_std_H palmDist_std_S{subID,3}];
    figure;
    b = bar(palmDist_avg);
    b(1).FaceColor = [22/255 160/255 133/255];
    b(2).FaceColor = [231/255 76/255 60/255];
    hold on
    
    % Finding the number of groups and the number of bars in each group
    ngroups = size(palmDist_avg,1);
    nbars = size(palmDist_avg,2);
    
    % Calculating the width for each bar group
    groupwidth = min(0.8,nbars/(nbars + 1.5));
    
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    % URL: https://uk.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x,palmDist_avg(:,i),palmDist_std(:,i),'k.');   
    end
    hold off
    xlabel('Session number');
    ylabel('Distance (mm)');
    legend(['Average Healthy'],['Subject ',num2str(subID)]);
    title('Total distance spanned by Palm Centroid during CPRT');
end

%% Method A: CoP DATA Averaging over each trial in one session for each subject

% Average CoP radial displacement for each session of each subject
for subType = 1:2
    subCount = subCount_cell{subType};
    subSesh = subSesh_cell{subType};
    trials = trials_cell{subType};

    for subID = subCount
        for session = subSesh
            avg_x = zeros(size(copData_x_norm{subType,subID,session},1), 1);
            avg_y = zeros(size(copData_y_norm{subType,subID,session},1), 1);
            avg_disp = zeros(size(copData_disp_norm{subType,subID,session},1), 1);

            for itt = 1:size(avg_disp,1)
                for taskID = 1:trials
                    avg_x(itt) = avg_x(itt) + copData_x_norm{subType,subID,session,taskID}(itt,1);
                    avg_y(itt) = avg_y(itt) + copData_y_norm{subType,subID,session,taskID}(itt,1);
                    avg_disp(itt) = avg_disp(itt) + copData_disp_norm{subType,subID,session,taskID}(itt,1);
                end % taskID
            end
            copData_x_avg_a{subType,subID,session} = avg_x/trials;
            copData_y_avg_a{subType,subID,session} = avg_y/trials;
            copData_disp_avg_a{subType,subID,session} = avg_disp/trials;
        end %session
        clear avg_x avg_y avg_disp;
    end %subID
end %subType

% First sample standard deviation of the average CoP radial displacement ^
for subType = 1:2
    subCount = subCount_cell{subType};
    subSesh = subSesh_cell{subType};
    trials = trials_cell{subType};
    
    for subID = subCount
        for session = subSesh
            std_x = zeros(size(copData_x_avg_a{subType,subID,session},1),1);
            std_y = zeros(size(copData_y_avg_a{subType,subID,session},1),1);
            std_disp = zeros(size(copData_disp_avg_a{subType,subID,session},1),1);
            for itt = 1:size(std_disp,1)
                for taskID = 1:trials
                    std_x(itt) = std_x(itt) + (copData_x_norm{subType,subID,session,taskID}(itt,1)-copData_x_avg_a{subType,subID,session}(itt,1))^2;
                    std_y(itt) = std_y(itt) + (copData_y_norm{subType,subID,session,taskID}(itt,1)-copData_y_avg_a{subType,subID,session}(itt,1))^2;
                    std_disp(itt) = std_disp(itt) + (copData_disp_norm{subType,subID,session,taskID}(itt,1)-copData_disp_avg_a{subType,subID,session}(itt,1))^2;
                end %taskID
            end
            if trials > 1
                copData_x_std_a{subType,subID,session} = sqrt(std_x/(trials-1));
                copData_y_std_a{subType,subID,session} = sqrt(std_y/(trials-1));
                copData_disp_std_a{subType,subID,session} = sqrt(std_disp/(trials-1));
            else
                copData_x_std_a{subType,subID,session} = 0;
                copData_y_std_a{subType,subID,session} = 0;
                copData_disp_std_a{subType,subID,session} = 0;
            end
        end %session
        clear std_x std_y std_disp;
    end %subID
end %subType

% Graph 1: CoP Radial Displacement Vs Task Completion Percentage
subH = 4;
for subS = [7,13]
    figure;
    task_complt = ((1:size(copData_disp_avg_a{1,subH},1))/(size(copData_disp_avg_a{1,subH},1)/100))';
    hold on;
    plot(task_complt, copData_disp_avg_a{1,subH,1}, 'Color', [22/255 160/255 133/255], 'LineWidth', 2);                                          % H, Mean
    plot(task_complt, copData_disp_avg_a{1,subH,1}-copData_disp_std_a{1,subH,1}, 'LineStyle', '-.', 'Color', [39/255 174/255 96/255], 'LineWidth', 1);    % H, Minimum
    plot(task_complt, copData_disp_avg_a{1,subH,1}+copData_disp_std_a{1,subH,1}, 'LineStyle', '--', 'Color', [41/255 128/255 185/255], 'LineWidth', 1);   % H, Maximum
    plot(task_complt, copData_disp_avg_a{2,subS,1}, 'Color', [231/255 76/255 60/255], 'LineWidth', 2);                                           % S, Mean
    plot(task_complt, copData_disp_avg_a{2,subS,1}-copData_disp_std_a{2,subS,1}, 'LineStyle', '-.', 'Color', [230/255 126/255 34/255], 'LineWidth', 1);   % S, Minimum
    plot(task_complt, copData_disp_avg_a{2,subS,1}+copData_disp_std_a{2,subS,1}, 'LineStyle', '--', 'Color', [192/255 57/255 43/255], 'LineWidth', 1);    % S, Maximum
    hold off;
    xlabel('Task Completion Percentage (%)');
    ylabel('CoP Displacement from starting position (mm)');
    title(['Avg Disp of CoP for H0', num2str(subH),' and S', num2str(subS),' during CPRT']);
    legend(['H0', num2str(subH),' (Mean)'], ['H0', num2str(subH),' (Mean - 1 SD)'], ['H0', num2str(subH),' (Mean + 1 SD)'], ...
            ['S', num2str(subS),' (Mean)'],  ['S', num2str(subS),' (Mean - 1 SD)'],  ['S', num2str(subS),' (Mean + 1 SD)'], 'Location', 'eastoutside');

    % Graph 2: CoP Displacement Vs Task Completion Percentage (New vs Old prosthetic)
    figure;
    hold on;
    plot(task_complt, copData_disp_avg_a{2,subS,1}, 'Color', [22/255 160/255 133/255], 'LineWidth', 2);                                          % H, Mean
    plot(task_complt, copData_disp_avg_a{2,subS,1}-copData_disp_std_a{2,subS,1}, 'LineStyle', '-.', 'Color', [39/255 174/255 96/255], 'LineWidth', 1);    % H, Minimum
    plot(task_complt, copData_disp_avg_a{2,subS,1}+copData_disp_std_a{2,subS,1}, 'LineStyle', '--', 'Color', [41/255 128/255 185/255], 'LineWidth', 1);   % H, Maximum
    plot(task_complt, copData_disp_avg_a{2,subS,2}, 'Color', [231/255 76/255 60/255], 'LineWidth', 2);                                           % S, Mean
    plot(task_complt, copData_disp_avg_a{2,subS,2}-copData_disp_std_a{2,subS,2}, 'LineStyle', '-.', 'Color', [230/255 126/255 34/255], 'LineWidth', 1);   % S, Minimum
    plot(task_complt, copData_disp_avg_a{2,subS,2}+copData_disp_std_a{2,subS,2}, 'LineStyle', '--', 'Color', [192/255 57/255 43/255], 'LineWidth', 1);    % S, Maximum
    hold off;
    xlabel('Task Completion Percentage (%)');
    ylabel('CoP Displacement from starting position (mm)');
    title(['Avg Disp of CoP for S', num2str(subS),' using Own (OP) and Advanced (AP) Prosthesis during CPRT']);
    legend('OP (Mean)', 'OP (Mean - 1 SD)', 'OP (Mean + 1 SD)', ...
           'AP (Mean)', 'AP (Mean - 1 SD)', 'AP (Mean + 1 SD)', 'Location', 'eastoutside');

    % Graph 3: CoP Displacement Vs Task Completion Percentage (W/ Vs W/out training)
    figure;
    hold on;
    plot(task_complt, copData_disp_avg_a{2,subS,2}, 'Color', [22/255 160/255 133/255], 'LineWidth', 2);                                          % H, Mean
    plot(task_complt, copData_disp_avg_a{2,subS,2}-copData_disp_std_a{2,subS,2}, 'LineStyle', '-.', 'Color', [39/255 174/255 96/255], 'LineWidth', 1);    % H, Minimum
    plot(task_complt, copData_disp_avg_a{2,subS,2}+copData_disp_std_a{2,subS,2}, 'LineStyle', '--', 'Color', [41/255 128/255 185/255], 'LineWidth', 1);   % H, Maximum
    plot(task_complt, copData_disp_avg_a{2,subS,3}, 'Color', [231/255 76/255 60/255], 'LineWidth', 2);                                           % S, Mean
    plot(task_complt, copData_disp_avg_a{2,subS,3}-copData_disp_std_a{2,subS,3}, 'LineStyle', '-.', 'Color', [230/255 126/255 34/255], 'LineWidth', 1);   % S, Minimum
    plot(task_complt, copData_disp_avg_a{2,subS,3}+copData_disp_std_a{2,subS,3}, 'LineStyle', '--', 'Color', [192/255 57/255 43/255], 'LineWidth', 1);    % S, Maximum
    hold off;
    xlabel('Task Completion Percentage (%)');
    ylabel('CoP Displacement from starting position (mm)');
    title(['Avg Disp of CoP for S', num2str(subS),' before (BT) and after (AT) Adv Prosthesis training during CPRT']);
    legend('BT (Mean)', 'BT (Mean - 1 SD)', 'BT (Mean + 1 SD)', ...
           'AT (Mean)', 'AT (Mean - 1 SD)', 'AT (Mean + 1 SD)', 'Location', 'eastoutside');
end %subSesh

%% Method A: Palm Centroid Disp Averaging over each trial in one session for each subject

% Average CoP radial displacement for each session of each subject
for subType = 1:2
    subCount = subCount_cell{subType};
    subSesh = subSesh_cell{subType};
    trials = trials_cell{subType};

    for subID = subCount
        for session = subSesh
            avg_disp = zeros(size(palmCentroidDispPerSub_norm{subType,subID,session},1), 1);
            for itt = 1:size(avg_disp,1)
                for taskID = 1:trials
                    avg_disp(itt) = avg_disp(itt) + palmCentroidDispPerSub_norm{subType,subID,session,taskID}(itt,1);
                end % taskID
            end
            palmCentroidDispPerSub_avg_a{subType,subID,session} = avg_disp/trials;
        end %session
        clear avg_disp;
    end %subID
end %subType

% First sample standard deviation of the average CoP radial displacement ^
for subType = 1:2
    subCount = subCount_cell{subType};
    subSesh = subSesh_cell{subType};
    trials = trials_cell{subType};
    
    for subID = subCount
        for session = subSesh
            std_disp = zeros(size(palmCentroidDispPerSub_avg_a{subType,subID,session},1),1);
            for itt = 1:size(std_disp,1)
                for taskID = 1:trials
                    std_disp(itt) = std_disp(itt) + (palmCentroidDispPerSub_norm{subType,subID,session,taskID}(itt,1)-palmCentroidDispPerSub_avg_a{subType,subID,session}(itt,1))^2;
                end %taskID
            end
            if trials > 1
                palmCentroidDispPerSub_std_a{subType,subID,session} = sqrt(std_disp/(trials-1));
            else
                palmCentroidDispPerSub_std_a{subType,subID,session} = 0;
            end
        end %session
        clear std_disp;
    end %subID
end %subType

% Graph 1: CoP Radial Displacement Vs Task Completion Percentage
subH = 4;
for subS = [7,13]
    figure;
    task_complt = ((1:size(palmCentroidDispPerSub_avg_a{1,subH},1))/(size(palmCentroidDispPerSub_avg_a{1,subH},1)/100))';
    hold on;
    plot(task_complt, palmCentroidDispPerSub_avg_a{1,subH,1}, 'Color', [22/255 160/255 133/255], 'LineWidth', 2);                                          % H, Mean
    plot(task_complt, palmCentroidDispPerSub_avg_a{1,subH,1}-palmCentroidDispPerSub_std_a{1,subH,1}, 'LineStyle', '-.', 'Color', [39/255 174/255 96/255], 'LineWidth', 1);    % H, Minimum
    plot(task_complt, palmCentroidDispPerSub_avg_a{1,subH,1}+palmCentroidDispPerSub_std_a{1,subH,1}, 'LineStyle', '--', 'Color', [41/255 128/255 185/255], 'LineWidth', 1);   % H, Maximum
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,1}, 'Color', [231/255 76/255 60/255], 'LineWidth', 2);                                           % S, Mean
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,1}-palmCentroidDispPerSub_std_a{2,subS,1}, 'LineStyle', '-.', 'Color', [230/255 126/255 34/255], 'LineWidth', 1);   % S, Minimum
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,1}+palmCentroidDispPerSub_std_a{2,subS,1}, 'LineStyle', '--', 'Color', [192/255 57/255 43/255], 'LineWidth', 1);    % S, Maximum
    hold off;
    xlabel('Task Completion Percentage (%)');
    ylabel('Palm Centroid Disp from starting position (mm)');
    title(['Avg Disp of Palm Centroid for H0', num2str(subH),' and S', num2str(subS),' during CPRT']);
    legend(['H0', num2str(subH),' (Mean)'], ['H0', num2str(subH),' (Mean - 1 SD)'], ['H0', num2str(subH),' (Mean + 1 SD)'], ...
            ['S', num2str(subS),' (Mean)'],  ['S', num2str(subS),' (Mean - 1 SD)'],  ['S', num2str(subS),' (Mean + 1 SD)'], 'Location', 'eastoutside');

    % Graph 2: CoP Displacement Vs Task Completion Percentage (New vs Old prosthetic)
    figure;
    hold on;
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,1}, 'Color', [22/255 160/255 133/255], 'LineWidth', 2);                                          % H, Mean
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,1}-palmCentroidDispPerSub_std_a{2,subS,1}, 'LineStyle', '-.', 'Color', [39/255 174/255 96/255], 'LineWidth', 1);    % H, Minimum
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,1}+palmCentroidDispPerSub_std_a{2,subS,1}, 'LineStyle', '--', 'Color', [41/255 128/255 185/255], 'LineWidth', 1);   % H, Maximum
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,2}, 'Color', [231/255 76/255 60/255], 'LineWidth', 2);                                           % S, Mean
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,2}-palmCentroidDispPerSub_std_a{2,subS,2}, 'LineStyle', '-.', 'Color', [230/255 126/255 34/255], 'LineWidth', 1);   % S, Minimum
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,2}+palmCentroidDispPerSub_std_a{2,subS,2}, 'LineStyle', '--', 'Color', [192/255 57/255 43/255], 'LineWidth', 1);    % S, Maximum
    hold off;
    xlabel('Task Completion Percentage (%)');
    ylabel('Palm Centroid Disp from starting position (mm)');
    title(['Avg Disp of Palm Centroid for S', num2str(subS),' using Own (OP) and Advanced (AP) Prosthesis during CPRT']);
    legend('OP (Mean)', 'OP (Mean - 1 SD)', 'OP (Mean + 1 SD)', ...
           'AP (Mean)', 'AP (Mean - 1 SD)', 'AP (Mean + 1 SD)', 'Location', 'eastoutside');

    % Graph 3: CoP Displacement Vs Task Completion Percentage (W/ Vs W/out training)
    figure;
    hold on;
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,2}, 'Color', [22/255 160/255 133/255], 'LineWidth', 2);                                          % H, Mean
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,2}-palmCentroidDispPerSub_std_a{2,subS,2}, 'LineStyle', '-.', 'Color', [39/255 174/255 96/255], 'LineWidth', 1);    % H, Minimum
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,2}+palmCentroidDispPerSub_std_a{2,subS,2}, 'LineStyle', '--', 'Color', [41/255 128/255 185/255], 'LineWidth', 1);   % H, Maximum
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,3}, 'Color', [231/255 76/255 60/255], 'LineWidth', 2);                                           % S, Mean
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,3}-palmCentroidDispPerSub_std_a{2,subS,3}, 'LineStyle', '-.', 'Color', [230/255 126/255 34/255], 'LineWidth', 1);   % S, Minimum
    plot(task_complt, palmCentroidDispPerSub_avg_a{2,subS,3}+palmCentroidDispPerSub_std_a{2,subS,3}, 'LineStyle', '--', 'Color', [192/255 57/255 43/255], 'LineWidth', 1);    % S, Maximum
    hold off;
    xlabel('Task Completion Percentage (%)');
    ylabel('Palm Centroid Disp from starting position (mm)');
    title(['Avg Disp of Palm Centroid for S', num2str(subS),' before (BT) and after (AT) Adv Prosthesis training during CPRT']);
    legend('BT (Mean)', 'BT (Mean - 1 SD)', 'BT (Mean + 1 SD)', ...
           'AT (Mean)', 'AT (Mean - 1 SD)', 'AT (Mean + 1 SD)', 'Location', 'eastoutside');
end %subSesh