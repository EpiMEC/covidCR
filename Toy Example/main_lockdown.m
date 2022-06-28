numSim = 1000;
% con medidas
direc = './data_lockdown/';
nameFile1 = 'totalS_';   nameFile2 = 'totalE_';   nameFile3 = 'totalR_';
nameFile4 = 'totalO1_';  nameFile5 = 'totalO2_';  nameFile6 = 'totalO3_';
nameFile7 = 'totalU1_';  nameFile8 = 'totalU2_';  nameFile9 = 'totalU3_';
nameFile10 = 'totalH1_'; nameFile11 = 'totalH2_'; nameFile12 = 'totalH3_';
nameFile13 = 'totalD1_'; nameFile14 = 'totalD2_'; nameFile15 = 'totalD3_';
nameFile16 = 'accumD1_'; nameFile17 = 'accumD2_'; nameFile18 = 'accumD3_';

for sim = 1:numSim
    [totalS,totalE,totalR,totalO1,totalO2,totalO3,totalU1,...
        totalU2,totalU3,totalH1,totalH2,totalH3,totalD1,...
        totalD2,totalD3,accumD1,accumD2,accumD3] = mainCall();
    parsave(sim,totalS,direc,nameFile1);
    parsave(sim,totalE,direc,nameFile2);
    parsave(sim,totalR,direc,nameFile3);
    parsave(sim,totalO1,direc,nameFile4);
    parsave(sim,totalO2,direc,nameFile5);
    parsave(sim,totalO3,direc,nameFile6);
    parsave(sim,totalU1,direc,nameFile7);
    parsave(sim,totalU2,direc,nameFile8);
    parsave(sim,totalU3,direc,nameFile9);
    parsave(sim,totalH1,direc,nameFile10);
    parsave(sim,totalH2,direc,nameFile11);
    parsave(sim,totalH3,direc,nameFile12);
    parsave(sim,totalD1,direc,nameFile13);
    parsave(sim,totalD2,direc,nameFile14);
    parsave(sim,totalD3,direc,nameFile15);
    parsave(sim,accumD1,direc,nameFile16);
    parsave(sim,accumD2,direc,nameFile17);
    parsave(sim,accumD3,direc,nameFile18);
end

function [totalS,totalE,totalR,totalO1,totalO2,totalO3,totalU1,...
    totalU2,totalU3,totalH1,totalH2,totalH3,totalD1,...
    totalD2,totalD3,accumD1,accumD2,accumD3] = mainCall()
tic
% parameters
T = 140;
numCounties = 15;
numInitExposed = 10;
load('nearCounties.mat') % counties connectivity
load('minmaxBubble.mat') % bubble sizes
min_coworkers = 5; delta_coworkers = 10;
min_sporadic  = 1; delta_sporadic  = 15;

coefSpor = .5;
coefFam = 1;
p  = .21;
pD = .13;
pM = .07;
pDM = 0.04;
pO = .10; % % isolated O
coefMin = 0.5;%0.012;
acum = readProbFile('prob');
familyFILE = 'setupFamily.csv';

%
pUseMask = .3;
minDist  = 10;

% Data.
% daily data (S, E, R: totals; O, U, H, D: per age group)
totalS  = zeros(1,T,'uint32'); accumD1 = zeros(1,T+1,'uint32');
totalE  = zeros(1,T,'uint32'); accumD2 = zeros(1,T+1,'uint32');
totalR  = zeros(1,T,'uint32'); accumD3 = zeros(1,T+1,'uint32');
totalO1 = zeros(1,T,'uint32'); totalU1 = zeros(1,T,'uint32');
totalH1 = zeros(1,T,'uint32'); totalD1 = zeros(1,T,'uint32');
totalO2 = zeros(1,T,'uint32'); totalU2 = zeros(1,T,'uint32');
totalH2 = zeros(1,T,'uint32'); totalD2 = zeros(1,T,'uint32');
totalO3 = zeros(1,T,'uint32'); totalU3 = zeros(1,T,'uint32');
totalH3 = zeros(1,T,'uint32'); totalD3 = zeros(1,T,'uint32');
% info per county
dataMap     = zeros(numCounties,T,'uint32');
dataMapAcum = zeros(numCounties,T,'uint32');
dataMapH    = zeros(numCounties,T,'uint32');
dataMapDead = zeros(numCounties,T,'uint32');
dataMapU    = zeros(numCounties,T,'uint32');
% 2. Create household network
[numIDs,numFam,idToCounty,idToFamily,idToAge,familyToIDs,countyToIDs] = setupFamilies(familyFILE);
% keep track of...
familyAdded = false(numFam,1); % ... added families
notadded = true(numIDs,1);     % ... not added IDs
% edge connectivity per layer
edges = cell(numIDs,3);
% number of contacts per ID: uniform value on [min, min+delta-1]
dgrLyr2 = uint8(min_coworkers+randi(delta_coworkers,numIDs,1)-1);
dgrLyr3 = uint8(min_sporadic +randi(delta_sporadic, numIDs,1)-1);
degrees = [dgrLyr2 dgrLyr3];
% list of possible contacts (to save time)
poolPeople   = cell(numCounties,1); % lista de elegibles por cantón
for k = 1:numCounties
    poolPeople{k} = uint32(cat(2,countyToIDs{nearCounties(k,:)>0}));
end
% 3. Initial conditions
% states: 0-S 1-E 2-O 3-U 4-H 5-R 6-D
stateID     = zeros(numIDs,1,'uint8');
timeStateID = zeros(numIDs,1,'uint8');
% numInitExposed initial exposed IDs
newExposedIDs = randi(numIDs,[numInitExposed,1]);
stateID(newExposedIDs) = 1; % S->E
timeStateID(newExposedIDs) = 1; % on day 1
% 4.
for day = 1:T
    % update network: include contacts for new exposed people
    for j = 1:numel(newExposedIDs)
        currID     = newExposedIDs(j);  % current ID
        numFamily  = idToFamily(currID);
        if(~familyAdded(numFamily)) % check if family has been added before
            familyNodes = familyToIDs{numFamily};  % IDs for family
            nodesToAdd = familyNodes(notadded(familyNodes));
            notadded(nodesToAdd) = false;          % update added nodes
            allEdges = combnk(familyNodes,2);      % family edges
            familyAdded(numFamily)= true;          % mark family as added
            for m = 1:size(allEdges,1)             % update edges
                edges{allEdges(m,1),1}(end+1) = allEdges(m,2);
                edges{allEdges(m,2),1}(end+1) = allEdges(m,1);
            end
        end
        %for layer = 2:3
        layer = 2;
        newConnections = max(0,degrees(currID,layer-1)-numel(edges{currID,layer}));
        if(newConnections>0)
            % who is eligible
            population = poolPeople{idToCounty(currID)};
            elegiblesW = population(randi(numel(population),1,3*newConnections));
            elegiblesW = setdiff(elegiblesW, [currID edges{currID,1} edges{currID,2} edges{currID,3}]);
            newCoworkers = elegiblesW(randperm(numel(elegiblesW),min(newConnections,numel(elegiblesW))));
            % find preexisting nodes on the graph
            nodesToAdd = newCoworkers(notadded(newCoworkers));
            % update added
            notadded(nodesToAdd) = false;
            % add edges
            allEdges = newCoworkers;
            % update neighborhs
            edges{currID,layer}(end+1:end+numel(allEdges)) = allEdges;
            for m = 1:numel(allEdges)
                edges{allEdges(m),layer}(end+1) = currID;
            end
        end
        %end
    end
    layer = 3;
    vector = find(stateID==2 | stateID==3)';
    for j = 1:numIDs
        edges{j,layer} = [];
    end
    for currID = vector
        newConnections = degrees(currID,layer-1);
        population = poolPeople{idToCounty(currID)};
        elegiblesW = population(randi(numel(population),1,2*newConnections));
        teem = [currID edges{currID,1} edges{currID,2} edges{currID,3}];
        %teeem = [currID cat(2,edges{currID,:})];
        elegiblesW = setdiff(elegiblesW,teem);
        newCoworkers = elegiblesW(randperm(numel(elegiblesW),min(numel(elegiblesW),newConnections)));
        % find preexisting nodes on the graph
        nodesToAdd = newCoworkers(notadded(newCoworkers));
        % update added
        notadded(nodesToAdd) = false;
        % add edges
        allEdges = newCoworkers;
        % update neighborhs
        edges{currID,layer} = allEdges;
        for m = 1:numel(allEdges)
            edges{allEdges(m),layer}(end+1) = currID;
        end
    end
    % probability of infection: state 2, 3, 4 can infect their contacts
    newExposedIDs = [];   % restart list of new E
    StoE = unique(cat(2,edges{(stateID==2|stateID==3|stateID==4),1:3}));
    % set daily % mask and social distancing
    if(day==20) % strict measures on day 20
        pUseMask = .9;
        minDist  = .7;
        minmaxBubble = randi(5,[numCounties 2]);
        minmaxBubble = sort(minmaxBubble,2);
    end
    if(day==80) % release measures on day 80
        pUseMask = .7;
        minDist  = .3;
        minmaxBubble = randi(15,[numCounties 2]);
        minmaxBubble = sort(minmaxBubble,2);
    end
    useMask  = randsample(2,numIDs,true,[1-pUseMask pUseMask]);
    pUseDist = (minDist+10*rand(numCounties,1))/100;
    for currID = StoE
        if(stateID(currID)==0)  % check susceptibles
            numCounty  = idToCounty(currID);
            sizeBubble = minmaxBubble(numCounty,1)+randi(minmaxBubble(numCounty,2)-minmaxBubble(numCounty,1)+1,1)-1;
            chosenN1 = edges{currID,2};
            chosenN2 = edges{currID,3};
            if(numel(chosenN1)+numel(chosenN2) > sizeBubble)
                nN = randi(sizeBubble,1);
                nN1 = min(nN,length(chosenN1));
                nN2 = min(sizeBubble-nN1,length(chosenN2));
                while(nN1+nN2<sizeBubble)
                    if(nN1<length(chosenN1))
                        nN1 = nN1 + 1;
                    end
                    if(nN2<length(chosenN2))
                        nN2 = nN2 + 1;
                    end
                end
                chosenN1 = chosenN1(randsample(length(chosenN1),nN1));
                chosenN2 = chosenN2(randsample(length(chosenN2),nN2));
            end
            % distribution for people that use mask and apply social distancing
            useMask1 = useMask(chosenN1);
            useMask2 = useMask(chosenN2);
            [~, useDist1] = histc(rand(numel(chosenN1),1),[0 1-pUseDist(numCounty) 1]);
            [~, useDist2] = histc(rand(numel(chosenN2),1),[0 1-pUseDist(numCounty) 1]);
            % 0-S 1-E 2-O 3-U 4-H 5-R 6-D
            states = stateID(edges{currID,1});             % layer 1
            iu = states == 3;
            id = states == 2 | states == 4;
            p1 = sum(iu) + pO*sum(id);
            risk = (1-coefFam*p).^(coefMin*p1);
            states = stateID(chosenN1);                    % layer 2
            iu = states == 3;
            id = states == 2 | states == 4;
            for k = 1:numel(iu)
                if(iu(k)==1)
                    if(useMask1(k)==2 && useDist1(k)==2)
                        risk = risk*((1-pDM).^coefMin);
                    elseif(useMask1(k)==2 && useDist1(k)==1)
                        risk = risk*((1-pM).^coefMin);
                    elseif(useMask1(k)==1 && useDist1(k)==2)
                        risk = risk*((1-pD).^coefMin);
                    else
                        risk = risk*((1-p).^coefMin);
                    end
                end
            end
            for k = 1:numel(id)
                if(id(k)==1)
                    if(useMask1(k)==2 && useDist1(k)==2)
                        risk = risk*((1-pDM).^(coefMin*pO));
                    elseif(useMask1(k)==2 && useDist1(k)==1)
                        risk = risk*((1-pM).^(coefMin*pO));
                    elseif(useMask1(k)==1 && useDist1(k)==2)
                        risk = risk*((1-pD).^(coefMin*pO));
                    else
                        risk = risk*((1-p).^(coefMin*pO));
                    end
                end
            end
            states = stateID(chosenN2);                    % layer 3
            iu = states == 3;
            id = states == 2 | states == 4;
            for k = 1:numel(iu)
                if(iu(k)==1)
                    if(useMask2(k)==2 && useDist2(k)==2)
                        risk = risk*((1-coefSpor*pDM).^coefMin);
                    elseif(useMask2(k)==2 && useDist2(k)==1)
                        risk = risk*((1-coefSpor*pM).^coefMin);
                    elseif(useMask2(k)==1 && useDist2(k)==2)
                        risk = risk*((1-coefSpor*pD).^coefMin);
                    else
                        risk = risk*((1-coefSpor*p).^coefMin);
                    end
                end
            end
            for k = 1:numel(id)
                if(id(k)==1)
                    if(useMask2(k)==2 && useDist2(k)==2)
                        risk = risk*((1-coefSpor*pDM).^(coefMin*pO));
                    elseif(useMask2(k)==2 && useDist2(k)==1)
                        risk = risk*((1-coefSpor*pM).^(coefMin*pO));
                    elseif(useMask2(k)==1 && useDist2(k)==2)
                        risk = risk*((1-coefSpor*pD).^(coefMin*pO));
                    else
                        risk = risk*((1-coefSpor*p).^(coefMin*pO));
                    end
                end
            end
            risk = 1 - risk; % probability infection
            [~, infected] = histc(rand,[0 risk 1]);
            if(infected==1)
                newExposedIDs = [newExposedIDs, currID]; %#ok<*AGROW>
            end
        end
    end
    % transitions
    timeStateID(newExposedIDs) = day; % first mark S->E
    stateID(newExposedIDs) = 1;
    elegible = find(stateID>0 & stateID<5); % all transitions but S->E
    states = stateID(elegible);             % all states
    indexD = [];                            % store new cases for accum data
    for j = 1:numel(states)
        state = states(j);
        ID = elegible(j);
        switch state
            case 1 % E
                time = day - timeStateID(ID);
                [~, event] = histc(rand,acum.E(time+1,:)); %#ok<*HISTC>
                if(event~=1)
                    stateID(ID) = event;
                    timeStateID(ID) = day;
                end
                if(event==2) % observed case
                    indexD = [indexD; ID];
                end
            case 2 % D
                time = day - timeStateID(ID);
                age  = idToAge(ID);
                switch age
                    case 1
                        [~, event] = histc(rand,acum.D1(time+1,:));
                    case 2
                        [~, event] = histc(rand,acum.D2(time+1,:));
                    case 3
                        [~, event] = histc(rand,acum.D3(time+1,:));
                end
                if(event~=2)
                    stateID(ID) = event;
                    timeStateID(ID) = day;
                end
            case 3 % U
                time = day - timeStateID(ID);
                [~, event] = histc(rand,acum.U(time+1,:));
                if(event~=3)
                    stateID(ID) = event;
                    timeStateID(ID) = day;
                end
            case 4 % H
                time = day - timeStateID(ID);
                age  = idToAge(ID);
                switch age
                    case 1
                        [~, event] = histc(rand,acum.HE1(time+1,:));
                    case 2
                        [~, event] = histc(rand,acum.HE2(time+1,:));
                    case 3
                        [~, event] = histc(rand,acum.HE3(time+1,:));
                end
                if(event~=4)
                    stateID(ID) = event;
                    timeStateID(ID) = day;
                end
        end
    end
    % save data
    is = (stateID==0)'; ie = (stateID==1)'; io = (stateID==2)';
    iu = (stateID==3)'; ih = (stateID==4)'; ir = (stateID==5)';
    id = (stateID==6)';
    accumD1(day+1)   = accumD1(day) + sum(idToAge(indexD)==1);
    accumD2(day+1)   = accumD2(day) + sum(idToAge(indexD)==2);
    accumD3(day+1)   = accumD3(day) + sum(idToAge(indexD)==3);
    totalS(day) = sum(is);  totalE(day) = sum(ie);  totalR(day) = sum(ir);
    % info according age groups
    ii = (idToAge==1);
    totalU1(day) = sum(iu & ii); totalO1(day) = sum(io & ii);
    totalH1(day) = sum(ih & ii); totalD1(day) = sum(id & ii);
    ii = (idToAge==2);
    totalU2(day) = sum(iu & ii); totalO2(day) = sum(io & ii);
    totalH2(day) = sum(ih & ii); totalD2(day) = sum(id & ii);
    ii = (idToAge==3);
    totalU3(day) = sum(iu & ii); totalO3(day) = sum(io & ii);
    totalH3(day) = sum(ih & ii); totalD3(day) = sum(id & ii);
    % save info per county
    acumCounty           = hist(idToCounty(indexD),1:numCounties)';
    dataMap(:,day+1)     = hist(idToCounty(io),1:numCounties); %#ok<*HIST,*UNRCH>
    dataMapAcum(:,day+1) = dataMapAcum(:,day)+uint32(acumCounty);
    dataMapH(:,day+1)    = hist(idToCounty(ih),1:numCounties);
    dataMapDead(:,day+1) = hist(idToCounty(id),1:numCounties);
    dataMapU(:,day+1)    = hist(idToCounty(iu),1:numCounties);
    %day
end
toc
end

function parsave(k,x,dir,nameFile)
nameFile = [dir '/' nameFile num2str(k) '.mat'];
save(nameFile, 'x')
end