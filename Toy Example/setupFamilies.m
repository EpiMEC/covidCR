function [numIDs,numFam,idToCounty,idToFamily,idToAge,...
    familyToIDs,countyToIDs] = setupFamilies(fileName)
% read demographic data (located in the same directory)
data = readtable(fileName);
houses_county  = data.Households;  age_group_01 = data.AG1;
people_county  = data.Population'; age_group_02 = data.AG2;
numCounties  = size(data,1);           % number of counties
numFam       = sum(houses_county); % number of households
% allocate arrays
people_families = zeros(1,numFam,'uint8'); % people per family
people_age      = cell(1,numCounties);     % age for ID per county
% create families based on Poisson distribution
index = 0;
for j = 1:numCounties
    local_houses     = houses_county(j);
    local_population = people_county(j);
    avrg_per_house   = local_population/local_houses;
    % create distribution for number of individuals per family
    number_people_per_house = poissrnd(avrg_per_house-1,local_houses,1)+1;
    % store number of people per family
    adv = numel(number_people_per_house);
    people_families(index+1:index+adv) = number_people_per_house;
    index = index + adv;
    % correct number of people per county (due to the random distribution)
    people_county(j) = sum(number_people_per_house);
    % include age groups
    n1 = floor(age_group_01(j)*people_county(j));
    n2 = floor(age_group_02(j)*people_county(j));
    n3 = people_county(j)-n1-n2;
    % create number of people per age labeled as 1, 2, 3
    age_groups = [ones(1,n1), 2*ones(1,n2), 3*ones(1,n3)];
    % sort age groups
    people_age{j} = age_groups(randperm(people_county(j)));
end
% accumulate sums for families and counties
cumsum_people_per_family = [0 cumsum(uint32(people_families))];
cumsum_people_per_county = [0 cumsum(people_county)];
numIDs = cumsum_people_per_county(end); % total population
% allocate arrays
idToAge     = zeros(1,cumsum_people_per_family(end),'uint8');
idToCounty  = zeros(1,cumsum_people_per_county(end),'uint8');
countyToIDs = cell(1,numCounties);
idToFamily  = zeros(1,cumsum_people_per_family(end),'uint32');
familyToIDs = cell(1,numFam);
% fill idToCounty, idToAge, countyToIDs
index = 0;
for j = 1:numCounties
    idToCounty(cumsum_people_per_county(j)+1:cumsum_people_per_county(j+1)) = j;
    countyToIDs{j} = cumsum_people_per_county(j)+1:cumsum_people_per_county(j+1);
    adv = numel(people_age{j});
    idToAge(index+1:index+adv) = people_age{j};
    index = index + adv;
end
% fill idToFamily, familyToIDs
for j = 1:numFam
    idToFamily(cumsum_people_per_family(j)+1:cumsum_people_per_family(j+1)) = j;
    familyToIDs{j} = uint32(cumsum_people_per_family(j)+1:cumsum_people_per_family(j+1));
end
end