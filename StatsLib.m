classdef StatsLib
    %StatsLib Functions to speed up stats hw and to learn
    %   Detailed explanation goes here
    
    % TODO:
    %   - Look through all functions and add annotations about input and
    %     output data formats
    %   - Make function titles consistent in style
	%		- GenerateObjectToGenerateHere_FROM_inputVar1ANDinputVar2...()
    %   - Make variable naming consistent
    %   - Change all internal function references to class.function() style
    %     so that it actually works because it is (Static)
    %   - Add new functions for what we learned:
    %       - Chebyshevsky
    %       - Box and Whisker
    %   - Let functions ingest each others data
    %   - Functions need ability to ingest continuous data not just integer data (done via sig fig mirroring for rounding)
    %   - Add variable input length via argc input ability and auto assume
    %   - Data ingest stream: array for acceptable data types of function data input
    %     Compare against the inputs and function outputs error if error present
    
    properties
        %Property1
    end
    
    methods
        function obj = untitled(inputArg1,inputArg2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    methods (Static)
        %Function Style Documentation:
        %   - [] around the output, for example like this: [output] indicates output is a vector
        %   - output formatted without [] like this: output = functionName(input) means output is a singular value


        %Helper Functions
            function isNonRepeating = CheckForRepeats(inputVect)
                %CheckForRepeats() self explanatory
                %	Inputs:
                %       - inputVect{vector<int>} => 1 by n array of int or floats
                %	Outputs:
                %       - isNonRepeating{bool}   => single bool output
                %	short explanation here
                isNonRepeating = false;

                for i=1:numel(inputVect)
                    if (isempty(find(inputVect(1:i) == inputVect(i))) == false)
                        break;
                    end
                    if i == numel(inputVect)
                        isNonRepeating = true;
                    end
                end
            end

            function [output] = GenerateRandomNonRepeatingIntegerArray(sizeArray)
                %GenerateRandomNonReatingIntegerArray() self explanatory
                %	Inputs:
                %       - sizeArray {int}      => 1 integer defining output array size
                %	Outputs:
                %       - output {vector<int>} => 1 by n integer array
                %	detailed explanation here
                
                output = zeros(1,sizeArray);

                for i = 1:sizeArray
                    generatedNum = randi(sizeArray);
                    while (isempty(find(output == generatedNum)) == false)
                        generatedNum = randi(sizeArray);
                    end

                    output(i) = generatedNum;
                end
            end
            
            function [bool] = CheckUserStringInputValid(inputValue, acceptedValuesArray)
                bool = false;
                for i = 1:numel(acceptedValuesArray)
                    %disp(acceptedValuesArray(i));
                    %disp(strcmp(inputValue,acceptedValuesArray(i)));
                    if strcmp(inputValue, acceptedValuesArray(i))
                        bool = true;
                        break;
                    end
                end
            end

            function [outputM] = TurnVectorIntoMatrix(dataVect, delineationCounter)
                %TurnVectorIntoMatrix() Generate matrix from input vector
                %	Takes input vector and creates new rows of equal length by adding proceeding values
                %   to list of items to cycle through
                %   row length = delineationCounter int value
                %	Inputs:
                %		- [dataVect]         {double, float, int} => n by 1 array
                %       - delineationCounter {int}                => int value indicating array row end point
                %	Outputs:
                %		-[outputM] {inputs.dataVect data type} => n/delineationCounter by delineationCounter matrix
                
                sizeM = [numel(dataVect)./delineationCounter, delineationCounter];
                outputM = zeros(sizeM);
                
                for i = 1:sizeM(1)
                    for j = 1:sizeM(2)
                        outputM(i,j) = dataVect((i-1).*delineationCounter + j);
                    end
                end
            end
            
            function [reorientedContainer] = CorrectOrientation(inputContainer, acceptedOrientations)
                %CorrectOrientation() | Checks input container orientation and fixes it
                %	Checks dimensions of container and reorients matrix or vector based on accepted value
                %	Inputs:
                %		-inputContainer{dataType} => 
                %	Outputs:
                %		-[reorientedContainer]{dataType} => 
                
                sizeCont = size(inputContainer);
                orientationCont = 0;
                % orientation map:
                %   Vect or Matrix of dimensions: R x C
                %   0 -> unresolved
                %   1 -> horizontal R < C
                %   2 -> vertical R > C
                %   3 -> square R = C

                if sizeCont(1) == sizeCont(2)
                    orientationCont = 3;
                else
                    if sizeCont(1) > sizeCont(2)
                        orientationCont = 1;
                    else
                        orientationCont = 2;
                    end
                end
                
                reorientedContainer = inputContainer;
                if numel(find(acceptedOrientations==orientationCont)) == 0
                    if (orientationCont == 1 || orientationCont == 2) && (acceptedOrientations == 3) 
                        % Turn into square by padding with zeros
                        newCont = zeros(max(sizeCont));
                        for i = 1:sizeCont(1)
                            for j = 1:sizeCont(2)
                                newCont(i,j) = inputContainer(i,j);
                            end
                        end
                        reorientedContainer = newCont;
                    else % Rotate
                        reorientedContainer = inputContainer';
                    end
                end
            end

        %Basic Data Manipulations
            function [outputList] = SortList_ascending(inputList)
                %SortList_ascending() sort left to right, least to most
                %   Inputs:
                %   Outputs:
                %   detailed description here
                outputList = sort(inputList,'ascend');
            end

            function [outputList] = SortList_descending(inputList)
                %SortList_descending() sort left to right, most to least
                outputList = sort(inputList,'descend');
            end
            
            function [outputList] = SortList(inputList, direction)
                %SortList() sort list based
                %	order list based on inputs
                %	Inputs:
                %		- inputList{vector<int, float, double>}  => number list or data array
                %	Outputs:
                %		- output{dataType} => 
                
                
                outputList = sort(inputList, direction);
            end

            function [output] = RandomizeArray(input)
                %RandomizeArray short explanation here
                %	Inputs:
                %	Outputs:
                %	detailed explanation here
                output = zeros(size(input));
                
                %generate random index array to reassign
                    randomIndices = StatsLib.GenerateRandomNonRepeatingIntegerArray(numel(input));

                for i = 1:numel(input)
                    output(i) = input(randomIndices(i));
                end
            end

            function [dataCopy] = RemoveDataPointRepeats(data)
                %FunctionName() short explanation here
                %	detailed description here
                %	Inputs:
                %		- input{dataType}  => 
                %	Outputs:
                %		- output{dataType} => 
                
                dataCopy = data;
                
                for i = 1:numel(data)
                    currentNumIndices = find(dataCopy == dataCopy(i));
                    
                    if numel(currentNumIndices) > 1
                        for j = 2:numel(currentNumIndices)
                            dataCopy(currentNumIndices(j)) = [];
                        end
                    end
                end
            end

        %Central Tendencies
            function avg = CalcMean(input)
                %CalcMean self explanatory
                %   Input:
                %       - input {int or vector<int>} => 
                %   Output:
                %       - avg {int} => 
                %   
                %   detailed explanation here
                avg = mean(input);
            end

            function avg = CalcAverage(input)
                avg = mean(input)
            end

            function avg = CalcAverage_inputList(input1,input2,varargin)
                allValsVect = zeros(1,numel(varargin)+2);
                allValsVect(1:2) = [input1, input2];
                
                for i = 1:numel(varargin)
                    allValsVect(i+2) = varargin{i};
                end

                avg = mean(allValsVect);
            end

            function medianVal = CalcMedian(data)
                medianVal = median(data);
            end

            function modeVal = CalcMode(data)
                dataHitCounts = zeros(size(data));
                for i = 1:numel(data)
                    dataHitCounts(i) = numel(find(data == data(i)));
                end

                if max(dataHitCounts) > 1
                    modeVal = mode(data);
                else
                    modeVal = [];
                    % add argc that makes this display turned on by default but can disable when nested
                    disp('data has no mode');
                end
            end

            function range = CalcDataRange(data)
                range = max(data) - min(data);
            end

            function midRangeVal = CalcDataMidRange(data)
                midRangeVal = (max(data) + min(data)) ./ 2;
            end

            function [output] = CalcCentralTendencies_ALL(inputArray)
                %CalcCentralTendencies() calculate most common central tendencies
                output = [mean(inputArray) median(inputArray) StatsLib.CalcMode(inputArray)];
            end

        %Average functions
            function arrayAvg = CalcAverage_Vect(data)
                %CalcDataAverage() calculate average of input array
                arrayAvg = mean(data);
            end
            
            function [output] = CalcAverages_MatrixRows(inputM)
                %FunctionName short explanation here
                %	Inputs:
                %	Outputs:
                %	detailed explanation here
                
                sizeInput = size(inputM);
                output = zeros(sizeInput(1),1);
                
                for i = 1:sizeInput(1)
                    output(i) = mean(inputM(i,:));
                end
            end

            function [output] = CalcAverage_Matrix(inputM)
                %FunctionName short explanation here
                %	Inputs:
                %	Outputs:
                %	detailed explanation here
                
                output = mean(StatsLib.CalcAverages_MatrixRows(inputM));
            end
            
            function weightedAvg = CalcWeightedMean_CategoryMeanVectAndWeightVect(catMeans, weightArray)
                %FunctionName() short explanation here
				%	detailed description here
				%	Inputs:
				%		- catMeans    {vector<float>} => 
				%		- weightArray {vector<float>} => 
				%	Outputs:
				%		- output{dataType} => 
				
				
				products = catMeans .* weightArray;
                weightedAvg = sum(products) ./ sum(weightArray);
            end
            
            function weightedAvg = CalcWeightedMean_AvgTableM(categoricalAvgTableM)
                weightedAvg = StatsLib.CalcWeightedMean(categoricalAvgTableM(:,1), categoricalAvgTableM(:,2));
            end
            
            function weightedAvg = CalcWeightedMean_CategoryDataMAndWeightVect(dataM, weightVect)
                sizeData = size(dataM);
                avgVect = zeros(sizeData(1),1);
                for i = 1:numel(avgVect)
                    avgVect(i) = mean(dataM(i,:));
                end
                
                weightedAvg = StatsLib.CalcWeightedMean(avgVect, weightVect);
            end
            
        %Measures of Variation - stdDev and Variance
            function output = CalcVariance_Population(data)
                %CalcVariance_Population() calculate population data variance
                %	detailed description here
                %	Inputs:
                %		-data{dataType} => 
                %	Outputs:
                %		-output{dataType} => 
                
                output = sum((data-mean(data)).^2) ./ numel(data);
            end

            function output = CalcVariance_Sample(data)
                %CalcVariance_Sample() calculate sample data variance
                %	detailed description here
                %	Inputs:
                %		-data{dataType} => 
                %	Outputs:
                %		-[output]{dataType} => 
                
                output = sum((data - mean(data)).^2) ./ (numel(data) - 1);
            end

            function output = CalcVariance(data, sampleType)
                %CalcVariance() calculate variance of either sample or population
                %	detailed description here
                %	Inputs:
                %		-data, sampleType{dataType} => 
                %	Outputs:
                %		-output{dataType} => 
                
                if sampleType == lower("population")
                    output = StatsLib.CalcVariance_Population(data);
                elseif sampleType == lower("sample")
                    output = StatsLib.CalcVariance_Sample(data);
                end
            end
            
            function stdDev = CalcStdDev_Population(data)
                %CalcStdDev_Population() short explanation here
                %	detailed description here
                %	Inputs:
                %		-data{dataType} => 
                %	Outputs:
                %		-stdDev{dataType} => 
                
                stdDev = sqrt(StatsLib.CalcVariance_Population(data));
            end

            function stdDev = CalcStdDev_Sample(data)
                %CalcStdDev_Sample() short explanation here
                %	detailed description here
                %	Inputs:
                %		-data{dataType} => 
                %	Outputs:
                %		-stdDev{dataType} => 
                
                stdDev = sqrt(StatsLib.CalcVariance_Sample(data));
            end

            function stdDev = CalcStdDev(data, sampleType)
                %CalcStdDev() short explanation here
                %	detailed description here
                %	Inputs:
                %		- data{dataType} => 
                %	Outputs:
                %       - stdDev{dataType} => 
                %       - sampleType
                
                if lower(sampleType) == "population"
                    stdDev = StatsLib.CalcStdDev_Population(data);
                elseif lower(sampleType) == "sample"
                    stdDev = StatsLib.CalcStdDev_Sample(data);
                end
            end

        %All measures of variation in a single function
            function [output] = CalcBaseDataVariationMeasures(data, sampleType) %modify to have default value for sampleType using nargin and varargin
                output = zeros(1,4);
                output(2) = mean(data);
                output(1) = range(data);
                if strcmp(sampleType,"population")
                    output(3) = StatsLib.CalcVariance_Population(data);
                    output(4) = sqrt(output(3));
                elseif strcmp(sampleType,"sample")
                    output(3) = StatsLib.CalcVariance_Sample(data);
                    output(4) = sqrt(output(3));
                else
                    error('Sample Type not accepted input value')
                end
            end

        %data binning / frequency distribution table functions
            %General Value Extraction
                function classWidth = CalcClassWidth(inputData, numClasses)
                    %TITLE short description
                    %   detailed explanation here
                    classWidth = ceil((max(inputData) - min(inputData))./ numClasses);
                end
                
                function [midpoints] = CalcClassMidpoints(MclassBounds)
                    midpoints = (MclassBounds(:,1) + MclassBounds(:,2)) ./ 2;
                end
                
                function tableMean = CalcMeanOfFreqTable_ClassBoundsAndVectFreqs(MclassBounds, vectFreqs)
                    tableMean = (sum(CalcClassMidpoints(MclassBounds).*vectFreqs)) ./ sum(vectFreqs);
                end
                
				function tableMean = CalcMeanofFreqTable_FreqTableM(freqTableM)
					classBounds = [freqTableM(:,1), freqTableM(:,2)];
					tableFreqs = freqTableM(:,3)';
					tableMean = StatsLib.CalcMeanOfFreqTable(classBounds, tableFreqs);
				end

                function binEdges = CalcClassEdges(classBoundsM)
                    sizeClassM = size(classBoundsM);
                    binEdges = zeros(sizeClassM(1)+1,1);
                    
                    edgeDelta = (classBoundsM(2,1) - classBoundsM(1,2)) ./ 2;
                    binEdges(1) = (classBoundsM(1,1) - edgeDelta);
                    for i = 1:sizeClassM(1)
                        binEdges(i+1) = classBoundsM(i,2) + edgeDelta;
                    end
                end

                function [cumulFreqs] = CalcCumulFreqs_FROM_FreqArray(freqArray)
                    cumulFreqs = freqArray;
                    for i = 2:numel(FreqArray)
                        cumulFreqs(i) = freqArray(i) + cumulFreqs(i-1);
                    end
                end

                function [counts] = CountFreqs_StringInput(dataVect, members)
                    %CountFreqs_StringInput() | Counts occurences of member vars in dataVect
                    %	detailed description here
                    %	Inputs:
                    %		-dataVect, members{dataType} => 
                    %	Outputs:
                    %		-[counts]{dataType} => 
                    
                    counts = zeros(size(members));
                    for i = 1:numel(members)
                        counts(i) = numel(find(dataVect == members(i)));
                    end
                end

                function [dataM] = CalcFreqs_StringInput(dataVect, members)
                    %CalcFreqs_StringInput() | Generate Frequency Table Values from String Array dataVect
                    %	generate numerical data table if you want string members to be displayed alongside you need a cell output rather than a matrix
                    %	Inputs:
                    %		-dataVect, members{dataType} => 
                    %	Outputs:
                    %		-[dataM]{dataType} => 
                    
                    freqs = StatsLib.CountFreqs_StringInput(dataVect, members);
                    dataM = [];
                end

            %Rebuild Data Given Freq Distribution Table
                function [generatedData] = GenerateDataVect_FreqTableVects(M_classBounds, vect_Freqs)
                    generatedData = zeros(1,sum(vect_Freqs));
            
                    index = 1;
                    for i = 1:numel(vect_Freqs)
                        classDelta = M_classBounds(i,2) - M_classBounds(i,1);
                        for j = 1:vect_Freqs(i)
                            generatedData(index) = M_classBounds(i,1) + randi(classDelta);
                            index = index + 1;
                        end
                    end
                end
            
                function [generatedData] = GenerateDataVect_FreqTableMatrix(freqTableM)
                    generatedData = StatsLib.GenerateDataVect_FreqTableVects([freqTableM(:,1), freqTableM(:,2)], freqTableM(:,3));
                end
                
            %Givens: Class bounds, raw data
                function [freqs] = CalcFreq_DataWithClassBounds(data, classBoundsM)
                    sizeClass = size(classBoundsM);
                
                    freqs = zeros(sizeClass(1),1);
                
                    freqs(1) = numel(find(data<=classBoundsM(1,2)));
                    for i = 2:sizeClass(1)
                        freqs(i) = numel(find(data<=classBoundsM(i,2))) - sum(freqs);
                    end
                end
                
                function [freqs, relFreqs, cumulFreqs] = CalcFreqs_DataWithClassBounds(data, classBoundsM)
                    freqs = StatsLib.CalcFreq_DataWithClassBounds(data, classBoundsM);
                    relFreqs = freqs ./ numel(data);
                    cumulFreqs = zeros(numel(freqs),1);
                    cumulFreqs(1) = freqs(1);
                    for i = 2:numel(cumulFreqs)
                        cumulFreqs(i) = cumulFreqs(i-1) + freqs(i);
                    end
                end
                
                function [dataM] = GenerateFreqTableMatrix_DataWithBins(classBoundsM, data)
                    dataM = [classBoundsM, CalcFreqsFromDataWithClassBounds(classBoundsM, data)];
                end
                
                function [dataM] = GenerateCompleteFreqTableMatrix_DataWithBins(classBoundsM, data)
                    [freqs, relFreqs, cumulFreqs] = StatsLib.CalcFreqs_DataWithClassBounds(classBoundsM, data);
                    midpoints = StatsLib.CalcClassMidpoints(classBoundsM);
                    
                    dataM = [classBoundsM, midpoints, freqs, relFreqs, cumulFreqs];
                end
                
            %Givens: data, numClasses
                function [generatedClasses, freqs, relFreqs, cumulFreqs] = GenerateFreqTableDataVects_Data(inputData, numClasses)
                    sizeData = numel(inputData);
                
                    %Prealloc
                        generatedClasses = zeros(numClasses, 2);
                        freqs = zeros(1, numClasses);
                        relFreqs = zeros(1,numClasses);
                        cumulFreqs = zeros(1,numClasses);
                    
                    classWidth = StatsLib.CalcClassWidth(inputData, numClasses);
                    
                    %Generate classes for binning
                        generatedClasses(1,:) = [min(inputData), min(inputData) + classWidth - 1];
                        for i = 2:numClasses
                            generatedClasses(i,:) = [generatedClasses(i-1,1) + classWidth, generatedClasses(i-1,1) + 2.*classWidth - 1];
                        end
                
                    %frequency Calcs
                        freqs(1) = numel(find(inputData<=generatedClasses(1,2)));
                        relFreqs(1) = freqs(1) ./ sizeData;
                        cumulFreqs(1) = freqs(1);
                        for i = 2:numClasses
                            freqs(i) = numel(find(inputData<=generatedClasses(i,2))) - cumulFreqs(i-1);
                            relFreqs(i) = freqs(i) ./ sizeData;
                            cumulFreqs(i) = cumulFreqs(i-1) + freqs(i);
                        end
                end
            
                function [dataCont] = GenerateFreqTableDataMatrix_Data(inputData, numClasses)
                    %class width arrya should be consolidated into another
                    %function or added as an extra addition to this
                    %functions inputs as a binary option to add it and what
                    %row to add it in, make this code more robust
                    
                    %classWidth = StatsLib.CalcClassWidth(inputData, numClasses);
                    %classWidthArray = zeros(numClasses, 1);
                    %for i = 1:numClasses
                    %    classWidthArray(i) = classWidth;
                    %end
                
                    [classes, freqs, relFreqs, cumulFreqs] = StatsLib.GenerateFreqTableDataVects_Data(inputData, numClasses);
                
                    dataCont = classes;
                    dataCont(:,3) = freqs;
                    dataCont(:,4) = relFreqs;
                    dataCont(:,5) = cumulFreqs;
                    %dataCont(:,6) = classWidthArray;
                end
            
                function [dataM] = GenerateFreqTableDataMatrix_AllDataVects(classBoundsM, freqs, relFreqs, cumulFreqs)
                    dataM = [classBoundsM, freqs, relFreqs, cumulFreqs];
                end
                
                function [dataM] = GenerateCompleteFreqTableDataMatrix_DataAndClassCount(inputData, numClasses)
                    [classes, freqs, relFreqs, cumulFreqs] = StatsLib.GenerateFreqTableDataVects_Data(inputData, numClasses);
                    midpoints = StatsLib.CalcClassMidpoints(classes);
                    dataM = (classes); 
                    dataM(:,3) = midpoints;
                    dataM(:,4) = freqs; 
                    dataM(:,5) = relFreqs;
                    dataM(:,6) = cumulFreqs;
                end
                
        %Generate MATLAB TABLE from data
			% // TODO: fix this shit, number of columns doesnt match input size, need auto recognition of size
            function [histogramTable] = GenerateFreqTableObject_FreqTableMatrix(freqTableDataM)
                %prep data with labels by putting it into arrays
                    Classes                = freqTableDataM(:,1:2);
                    Midpoints              = freqTableDataM(:,3);
                    Frequencies            = freqTableDataM(:,4);
                    Relative_Frequencies   = freqTableDataM(:,5);
                    Cumulative_Frequencies = freqTableDataM(:,6);
                    Class_Widths           = freqTableDataM(:,7);
                
                %generate table
                    histogramTable = table(...
                        Classes,...
                        Midpoints,...
                        Frequencies,...
                        Relative_Frequencies,...
                        Cumulative_Frequencies,...
                        Class_Widths...
                    );
            end
            
            function [histogramTable] = GenerateFreqTableObject_DataAndClassCount(data, numClasses)
                dataM = GenerateCompleteFreqTableDataMatrix_Data(data, numClasses);
                histogramTable = GenerateCompleteFreqTable_DataMatrix(dataM);
            end
            
        %Histogram plot function here
            %Matlab not binning, prebinned via defined classes and provided
            %frequencies
                %class bounds & freqs
                    function [h] = GenerateHistogram_BinsFreqs(classBoundsM, freqs)
                        classEdges = CalcClassEdges(classBoundsM);
                        h = histogram('BinEdges',classEdges,'BinCounts',freqs);
                    end
            
                %freq dist table data
                    function [h] = GenerateHistogram_FreqTableMatrix(freqTableM)
                        binCounts = [];
                        midPoints = StatsLib.CalcClassMidpoints([freqTableM(:,1) freqTableM(:,2)]);
                        
                        %figure out normal midpoint placement here and
                        %standardize the placement across the rest of the
                        %functions
                        if midPoints == freqTableM(:,3) %normal
                            binCounts = freqTableM(:,4);
                        elseif midPoints == freqTableM(:,4)
                            binCounts = freqTableM(:,3);
                        else
                            binCounts = freqTableM(:,3);
                        end
                        %conventionally midpoints should be found in table
                        %column 3 right after the class boundary LV and UV
                        %values
                        
                        h = histogram('BinEdges',CalcClassEdges([freqTableM(:,1), freqTableM(:,2)]), 'BinCounts', binCounts);
                    end
            
            %Matlab data processing happens b/c raw data input
                %data , bins
                    function [h] = GenerateHistogram_DataAndBins(data, classBoundsM)
                        h = histogram(data, CalcClassEdges(classBoundsM));
                    end
        
        %Measures of Position - Quartiles and Outliers
            %Percentile Functions
                function percentileVal = CalcPercentileOfValue(val, data)
                    %CalcPercentileOfValue() | Takes dataset and value and calculates percentile value corresponding to input data value
                    %	detailed description here
                    %	Inputs:
                    %		-val, data{dataType} => 
                    %	Outputs:
                    %		-percentileVal{dataType} => 
                    
                    dataSorted = sort(data);
                    percentileVal = (numel(find(dataSorted<val)) ./ numel(data)) .* 100;
                end

                function percentileIndex = CalcPercentileIndex(percentileWeWant, data)
                    %CalcPercentileIndex() | Finds index value of percentile value we are looking for within the data
                    %	detailed description here
                    %	Inputs:
                    %		-percentileWeWant, data{dataType} => 
                    %	Outputs:
                    %		-percentileIndex{dataType} => 
                    
                    if percentileWeWant > 0 && percentileWeWant < 1
                        percentileUsing = percentileWeWant .* 100;
                    else
                        percentileUsing = percentileWeWant;
                    end

                    percentileIndex = round((percentileUsing ./ 100) .* numel(data),0) + 1;
                end

                function valueCorresponding2Percentile = CalcValueOfPercentile(percentileWeWant, data)
                    %CalcPercentileValueFromData() | Outputs value in dataset corresponding to percentile value we want
                    %	detailed description here
                    %	Inputs:
                    %		-percentileWeWant, data{dataType} => 
                    %	Outputs:
                    %		-valueCorresponding2Percentile{dataType} => 
                    
                    dataSorted = sort(data);
                    valueCorresponding2Percentile = dataSorted(StatsLib.CalcPercentileIndex(percentileWeWant, data));
                end

            %Locator formula
                function medianIndex = CalcMedianIndex(dataset)
                    %CalcMedianIndex() calculates index of median from dataset
                    %	detailed description here
                    %	Inputs:
                    %		- input{dataType}  => 
                    %	Outputs:
                    %		- output{dataType} => 
                    
                    datasetLength = numel(dataset);
                    medianIndex = (datasetLength + 1) ./ 2;
                end

            %Quartile Functions
                function [output] = CalcQuartiles(data)
                    %CalcQuartiles() Calculate the quartile data values from the input data
                    %	Calculate min, q1, q2, q3, max values for quartile data ranges
                    %	Inputs:
                    %		-data{dataType} => 
                    %	Outputs:
                    %		-output{dataType} => 
                
                    data_sorted = sort(data);

                    dataMedian = median(data);
                    medianIndex = StatsLib.CalcMedianIndex(data);
                    
                    dataQ1 = -1;
                    dataQ3 = -1;

                    if rem(numel(data),2) == 0
                        dataQ1 = median(data_sorted(1:floor(medianIndex)));
                        dataQ3 = median(data_sorted(ceil(medianIndex):numel(data)));
                    else
                        dataQ1 = median(data_sorted(1:medianIndex));
                        dataQ3 = median(data_sorted(medianIndex:numel(data)));
                    end

                    %if rem(numel(data),2) == 0
                    %    %disp('even');
                    %    dataQ1 = median(data_sorted(1:(numel(data)./2)));
                    %    dataQ3= median(data_sorted((numel(data)./2+1):numel(data)));
                    %else
                    %    disp('odd');
                    %    %middle value is excised from data set
                    %    middleIndex = find(dataMedian==data_sorted);
                    %    dataQ1 = median(data_sorted(1:middleIndex-1));
                    %    dataQ3 = median(data_sorted(middleIndex+1:numel(data)));
                    %end

                    output = [min(data), dataQ1, dataMedian, dataQ3, max(data)];
                end
            
                function [output] = CalcIQR(data)
                    %CalcIQR() calculate the distance between q3 and q1 of a dataset
                    %	detailed description here
                    %	Inputs:
                    %		-data{dataType} => 
                    %	Outputs:
                    %		-output{dataType} => 
                
                    data_quartile = StatsLib.CalcQuartiles(data);
                    output = data_quartile(4) - data_quartile(2);
                end

                function [output] = CalcQuartileMetrics_ALL(data)
                    %CalcQuartileData_ALL() All relevant quartile data output as a matrix
                    %	Quartile data in addition to IQR and other quartile information
                    %	Inputs:
                    %		-data{dataType} => 
                    %	Outputs:
                    %		-output{dataType} => 
                
                    output = StatsLib.CalcQuartiles(data);
                    output(6) = output(4) - output(2);
                end
                
                function [output] = CalcQuartileMetrics_IQRReadableWithIQR(data)
                    container = StatsLib.CalcQuartileMetrics_ALL(data);
                    values = zeros(1,8);
                    values(1) = container(1); %min
                    values(2) = container(2) - container(6); % lower outlier threshold
                    values(3) = container(2); % Q1
                    values(4) = container(3); % median
                    values(5) = container(4); % Q3
                    values(6) = container(4) + container(6); % upper outlier threshold
                    values(7) = container(5); % max
                    values(8) = container(6); % IQR
                    %disp('ran');

					labels = {"min", "lower outlier threshold", "Q1", "Q2 (median)", "Q3", "upper outlier threshold", "max", "IQR"};
					%output(2,:) = values;
					disp(labels);
					output = values;
				end

            %Outlier Functions
                function [lowerOutliers, upperOutliers] = CalcOutliersViaIQR(data)
                    %CalcOutliers() short explanation here
                    %	detailed description here
                    %	Inputs:
                    %		-data{dataType} => 
                    %	Outputs:
                    %		-[outliers]{dataType} => 
                
                    data_metrics_quartile = StatsLib.CalcQuartileMetrics_ALL(data);
                
                    boundBuffer = data_metrics_quartile(6) .* 1.5;
                    bounds = [data_metrics_quartile(2)-boundBuffer, data_metrics_quartile(4)+boundBuffer];

                    lowerOutlierIndices = find(data < bounds(1));
                    upperOutlierIndices = find(data > bounds(2));

                    lowerOutliers = zeros(1,numel(lowerOutlierIndices));
                    for i = 1:numel(lowerOutlierIndices)
                        lowerOutliers(i) = data(lowerOutlierIndices(i));
                    end

                    upperOutliers = zeros(1,numel(upperOutlierIndices));
                    for i = 1:numel(upperOutlierIndices)
                        upperOutliers(i) = data(upperOutlierIndices(i));
                    end
                end

                function [outliers] = CalcOutliersViaIQR_Consolidated(data)
                    %CalcOutliers_Consolidated() short explanation here
                    %	detailed description here
                    %	Inputs:
                    %		-input{dataType} => 
                    %	Outputs:
                    %		-output{dataType} => 
                
                    [a, b] = StatsLib.CalcOutliersViaIQR(data);
                    outliers = [a, b];
                end
            
            %All values together in one function
                function [quartileInfo, outliers] = CalcQuartileDataAndOutlierViaIQRData(data)
                    quartileInfo = StatsLib.CalcQuartileMetrics_IQRReadableWithIQR(data);
                    outliers = StatsLib.CalcOutliers_Consolidated(data);
                end
            
        %Chebyshevsky Functions
            function [output] = CalcPercentDataInStdDevRangeViaChebyshevsky(n_StdDevs)
                %CalcChebyshevsky() short explanation here
                %	detailed description here
                %	Inputs:
                %		-input{dataType} => 
                %	Outputs:
                %		-[output]{dataType} => 
                
                output = 1 - (1 ./ (n_StdDevs).^2);
            end

            function output = CalcPercentDataInStdDevRangeViaChebyshevsky_RawData(rangeYouWant, sampleMean, sampleStdDev)
                %CalcPercentDataInStdDevRangeViaChebyshevksy_RawData() short explanation here
                %	detailed description here
                %	Inputs:
                %		- input{dataType}  => 
                %	Outputs:
                %		- output{dataType} => 
                
                deltas = rangeYouWant - sampleMean;
                %deltasInStdDevs = deltas ./ sampleStdDev;

                %percentVal = StatsLib.CalcPercentDataInStdDevRangeViaChebyshevsky(deltaInStdDevs);
                % ^ This does not work because the percent width is calculated indepedently vs actual width
                
                %Needs to be compensated for syhmmetric distribution even though chebyshevsky
                %compensates for kurtosis?

                k = 0;
                if deltas(1) == deltas(2)
                    k = deltas(1) ./ stdDev;
                end

                if k > 1
                    output = StatsLib.CalcPercentDataInStdDevRangeViaChebyshevsky(k);
                else
                    output = -1;
                end
            end

            function output = CalcNumDataInStdDevRangeViaChebyshevsky_RawData(rangeYouWant, sampleMean, sampleStdDev)
                %CalcNumDataInStdDevRangeViaChebyshevsky_RawData() short explanation here
                %	detailed description here
                %	Inputs:
                %		-rangeYouWant, sampleMean, sampleStdDev{dataType} => 
                %	Outputs:
                %		-output{dataType} => 
                
                
            end

        %Normal Distribution Functions
            function [empiricalRuleM] = OutputEmpiricalRuleInfoM()
                empiricalRuleM = {...
                    "-n" -3 -2 -1 0 1 2 3 "n";...
                    0.15 2.35 13.5 34 0 34 13.5 2.35 0.15};
            end

            function [normalDistBoundaries] = CalcNormalDistBorderValues(data, sampleSet) 
                %need to change dataType to some other variable name to describe sampling distribution widthjhnhjkjhkkjh
                %add normal distribution shit here
                dataMean = mean(data);
                dataStdDev = StatsLib.CalcStdDev(data, sampleSet);

                normalDistBoundaries = [-3 -2 -1 0 1 2 3];
                normalDistBoundaries = normalDistBoundaries .* dataStdDev;
                normalDistBoundaries(4) = dataMean;
            end

			function output = FindOutliersUsingEmpiricalRule_SamplesAndMuAndSigmaAndThreshold(data, mu, sigma, thresholdVal, varargin)
				output = [];
				lowerOutliers = [];
				upperOutliers = [];
				if nargin > 4
				else
					upperBound = mu + thresholdVal .* sigma;
					lowerBound = mu - thresholdVal .* sigma;
					upperOutlierIndices = find(data > upperBound);
					lowerOutlierIndices = find(data < lowerBound);

					%lowerOutliers = [1];
					%upperOutliers = [];
					lowerOutliers = data(lowerOutlierIndices);
					upperOutliers = data(upperOutlierIndices);
				end

				output = [lowerOutliers, upperOutliers];
			end

			function [output] = FindOutliersUsingEmpircalRule(data, varargin)
				
			end

        %Skewing Functions and Kurtosis
            function skewDirection = DetermineDataSkew(data)
                %DetermineDataSkew() evaluate skew or tail difference from symmetric distribution
                %	detailed description here
                %	Inputs:
                %		-data{dataType} => 
                %	Outputs:
                %		-skewDirection{dataType} => 
                
                dataMean = mean(data);
                dataMedian = median(data);

                if dataMean < dataMedian
                    skewDirection = -1;
                elseif dataMean > dataMedian
                    skewDirection = 1;
                else
                    skewDirection = 0;
                end
            end

            function centralTendencyName = DetermineBestCentalTendency(data)
                %DetermineBestCentalTendency() Evaluate best measure of central tendency based on data
                %	detailed description here
                %	Inputs:
                %		-data{dataType} => 
                %	Outputs:
                %		-centralTendencyName{dataType} => 
                
                centralTendencyName = "null";

                if isnumeric(data(randi(numel(data))))
                    skewDirection = StatsLib.DetermineDataSkew(data);
                    if skewDirection == 0
                        centralTendencyName = "mean";
                    else
                        centralTendencyName = "median";
                    end
                else
                    centralTendencyName = "mode";
                end
            end

            %kurtosis function here

        % C4 - Probability Distributions
            %Probability Histogram Functions
                function isValid = VerifyProbabilityDistribution_FROM_probVect(probVect)
                    isValid = false;
                    freqSum = 0;

                    for i = 1:numel(probVect)
                        if probVect(i) > 1 || probVect(i) < 0
                            break;
                        end

                        freqSum = freqSum + probVect(i);

                        if freqSum > 1
                            break;
                        end

                        if i == numel(freqSum)
                            isValid = true;
                        end
                    end
                end

                function isValid = VerifyProbabilityDistribution_FROM_probMatrix(probMatrix)
                    %VerifyProbabilityDistribution() short explanation here
                    %	detailed description here
                    %	Inputs:
                    %		-probMatrix{dataType} => 
                    %	Outputs:
                    %		-isValid{dataType} => 
                    
                    isValid = false;
                    freqSum = 0;

                    for i = 1:numel(probMatrix(1,:))
                        if probMatrix(2,i) > 1 || probMatrix(2,i) < 0
                            break;
                        end

                        freqSum = freqSum + probMatrix(2,i);

                        if freqSum > 1
                            break;
                        end

                        if i == numel(probMatrix(1,:))
                            isValid = true;
                        end
                    end
                end
            
                function [output] = GenerateProbDistTableM_FROM_OutcomeVectANDProbVect(vect_Outcomes, vect_Probs)
                    %GenerateProbabilityDistributionTable_Matrix() short explanation here
                    %	detailed description here
                    %	Inputs:
                    %		-input{dataType} => 
                    %	Outputs:
                    %		-[output]{dataType} => 
                    
                    if numel(vect_Outcomes) == numel(vect_Probs)
                        output = vect_Outcomes';
                        output(:,2) = vect_Probs';
                    else
                        output = [-1,-1];
                    end
                end

                function [output] = GenerateProbDistTableM_FROM_FreqTableMatrix(freqTableM)
                    %GenerateProbabilityDistributionTable_Matrix() short explanation here
                    %	detailed description here
                    %	Inputs:
                    %       needs to be updated to accept all the possible freqTableM
                    %		-input{dataType} => 
                    %	Outputs:
                    %		-[output]{dataType} => 
                    
					% TODO: DATA ENTRY FORMAT NEEDS REVIEW FORGET WHAT COLUMN IS WHAT IN TEH FREQUNECY TABLE INPUTS
                    output = StatsLib.GenerateProbDistTableMatrix_FROM_OutcomeVectANDProbVect(FreqTableM(:,1),(FreqTableM(:,2)./sum(FreqTableM(:,2))));
                end

                function [output] = GenerateHistogram_ProbDistTable(DataMatrix_ProbDistTable)
                    %GenerateHistogram_ProbabilityDistributionTable() short explanation here
                    %	detailed description here
                    %	Inputs:
                    %		-{dataType} => 
                    %	Outputs:
                    %		-[output]{dataType} => 
                    
                    %Main Code Here
                end

            %Probability Distribution Central Values
                function output = CalcProbDist_Mean(probTableM)
                    %Calc_ProbDistribution_Mean() short explanation here
                    %	detailed description here
                    %	Inputs:
                    %		-probTableM {vector<vector<float>>} => vertically aligned two row table
					%											   column 1: x Values i.e. possible outcomes
					%											   column 2: P(x) values
                    %	Outputs:
                    %		-output {float} => mean of probabilty table
                    sizeInput = size(probTableM);

					if sizeInput(1) > sizeInput(2)
                        output = sum(probTableM(:,1).*probTableM(:,2));
                    else
                        output = sum(probTablem(1,:).*probTableM(2,:));
                    end
                end

                function output = CalcProbDist_Variance(probTableM)
                    %Calc_ProbDistribution_Variance() short explanation here
                    %	detailed description here
                    %	Inputs:
                    %		-input{dataType} => 
                    %	Outputs:
                    %		-[output]{dataType} => 
                    
                    currentMean = StatsLib.CalcProbDist_Mean(probTableM);
                    output = sum((probTableM(:,1) - currentMean).^2 .* probTableM(:,2));
                end
                
                function output = CalcProbDist_StdDev(input)
                    %Calc_ProbDist_StdDe() short explanation here
                    %	detailed description here
                    %	Inputs:
                    %		-input{dataType} => 
                    %	Outputs:
                    %		-[output]{dataType} => 
                    
                    output = sqrt(StatsLib.CalcProbDist_Variance(input));
                end

                function stdDevFreqTable = CalcStdDev_OFFreqTable(classBounds, classFreqs)
                    %CalcProbDist_FreqTableStdDev() | short explanation here
                    %	detailed description here
                    %	Inputs:
                    %		-classBounds, classFreqs{dataType} => 
                    %	Outputs:
                    %		-stdDevFreqTable{float} => 
                
                    % THIS IS FOR SAMPLES NOT POPULATIONS

                    classMidpoints = StatsLib.CalcClassMidpoints(classBounds);
                    numSamples = sum(classFreqs);
                    sumFreqTimesMidpoint = sum(classMidpoints .* classFreqs);
                    sumFreqTimesMidpointSquared = sum(classMidpoints.^2 .* classFreqs);
                    varianceValue = (sumFreqTimesMidpointSquared - (sumFreqTimesMidpoint.^2 ./ numSamples)) ./ (numSamples - 1); %change this to numSamples from (numSamples - 1) for populations
                    stdDevFreqTable = sqrt(varianceValue);
                end

                function [output] = CalcProbDist_AllStatMeasures(input)
                    %Calc_ProbDist_AllStatMeasures() short explanation here
                    %	detailed description here
                    %	Inputs:
                    %		-input{dataType} => 
                    %	Outputs:
                    %		-[output]{dataType} => 
                    
                    output = zeros(1,3);
                    output(1) = StatsLib.CalcProbDist_Mean(input);
                    output(2) = StatsLib.CalcProbDist_Variance(input);
                    output(3) = StatsLib.CalcProbDist_StdDev(input);
                end
                
                function coeffOfVariation = CalcCoefficientOfVariation(data, varargin)
                    %CalcCoefficientOfVariation() | Calculates Coefficient of Variation
                    %	Coefficient of variation is a risk 2 reward metric for different data series
                    %   Formula: stdDev / mean x 100 - used to understand spread of data set
                    %	Inputs:
                    %		-data{vector<int,float,double>} => n x 1 or 1 x n | data containing discrete and/or continuous data
                    %       -varagin{
                    %                sampleSpace{string} => string containing either "population" or "sample" to declare sample space for corresponding formula to be used
                    %               }
                    %	Outputs:
                    %		-coeffOfVariation{int,float,double} => value of coeffiecient of variation [unit: percent]
                    
                    if nargin == 2
                        acceptedInputs = ["population", "sample"];
                        if StatsLib.CheckUserStringInputValid(lower(varargin{1}), acceptedInputs)
                            if strcmp(varargin{1},"population")
                                coeffOfVariation = StatsLib.CalcStdDev_Population(data) ./ mean(data) .* 100; % in percent
                            elseif strcmp(varargin{1},"sample")
                                coeffOfVariation = StatsLib.CalcStdDev_Sample(data) ./ mean(data) .* 100; % in percent
                            end
                        end
                    elseif nargin == 1 %Default to sample
                        coeffOfVariation = StatsLib.CalcStdDev_Sample(data) ./ StatsLib.CalcMean(data) .* 100; % in percent
                    else
                        error('Too many inputs!');
                    end
                end

        % 4.2 - Binomial Distributions
            %Calculate Binomial Probabilities
                function output = CalcBinomialProbability_base(p,n,x)
                    output = nchoosek(n,x) .* p.^x .* (1-p).^(n-x);
                end

                function [output] = CalcBinomialProbability(p,n,x, varargin)
                    %CalcBinomialProbability() short explanation here
                    %	detailed description here
                    %	Inputs:
                    %		-p,n,x{dataType} =>
                    %       -calcDirection {string} => number of probabilities to from x position,
                    %                                  can be: '=', '<', '>'
                    %	Outputs:
                    %		-[output]{dataType} => 

                    %Input Arg Parsing
                        validInputs = ["=", ">", "<", ">=", "<="];
                        if nargin == 4
                            probWidth = varargin{1};
                            if StatsLib.CheckUserStringInputValid(probWidth, validInputs) == false
                                error('Incorrect comparison arg(varargin{4})')
                            end
                        elseif nargin == 3
                            probWidth = '=';
                        elseif nargin > 4
                            error('Too many input arguments. Max 4 arguments.')
                        elseif nargin < 3
                            error('Not enough inputs. Minimum 3 input arguments.');
                        else
                            error('unknown error');
                        end

                    %Calcs
                        if strcmp(probWidth,'=')
                            output = StatsLib.CalcBinomialProbability_base(p,n,x);
                        elseif strcmp(probWidth,'<')
                            output = 0;
                            for i = 0:x-1
                                output = output + StatsLib.CalcBinomialProbability_base(p,n,i);
                            end
                        elseif strcmp(probWidth,'>')
                            output = 0;
                            for i = x+1:1:n
                                output = output + StatsLib.CalcBinomialProbability_base(p,n,i);
                            end
                        elseif strcmp(probWidth,'>=')
                            output = 0;
                            for i = x:1:n
                                output = output + StatsLib.CalcBinomialProbability_base(p,n,i);
                            end
                        elseif strcmp(probWidth,'<=')
                            output = 0;
                            for i = 0:x
                                output = output + StatsLib.CalcBinomialProbability_base(p,n,i);
                            end
                        end
                end
                
            %Binomial Probability between values
                function probVal = CalcBinomialProbabilityBetweenXvals(p,n,xVal1,xVal2)
                    probVal = 0;
                    xVals = [xVal1, xVal2];
                    for i = min(xVals):1:max(xVals)
                        probVal = probVal + StatsLib.CalcBinomialProbability(p,n,i);
                    end
                end

			%Generate Binomial Probability Table
				function [binomialTableM] = GenerateBinomialDistributionTable(p,n,xValueVect)
					%FunctionName() short explanation here
					%	detailed description here
					%	Inputs:
					%		- input{dataType}  => 
					%	Outputs:
					%		- output{dataType} => 
					
					binomialTableM = xValueVect';

					for i = 1:numel(xValueVect)
						binomialTableM(i,2) = StatsLib.CalcBinomialProbability(p,n,xValueVect(i));
					end
				end

            %Binomial Central tendencies
                function binomAvg = CalcBinomialMean(p,n)
                %CalcBinomialMean - calculate binomial distribution mean
                %
                % Syntax: binomAvg = CalcBinomialMean(p,n)
                %
                % Long description
                    binomAvg = n .* p;
                end

                function output = CalcBinomialVariance(p,n)
                    output = n .* p .* (1-p);
                end

                function output = CalcBinomialStdDev(p,n)
                    output = sqrt(StatsLib.CalcBinomialVariance(p,n));
                end

                function [output] = CalcBinomialCentralTendencies_ALL(p,n)
                    output = zeros(1,3);
                    output(1) = StatsLib.CalcBinomialMean(p,n);
                    output(2) = StatsLib.CalcBinomialVariance(p,n);
                    output(3) = StatsLib.CalcBinomialStdDev(p,n);
                end

        % C5 - Normal Distributions
            % 5.1 - Gaussain or Normal Distributions
                %Calc Z-score
                    function zScore = CalcZScore(x,mu,stdDev)
                        zScore = (x - mu) ./ stdDev;
                    end
                    
                %Compare Z-Scores
                    function [zScoresVect] = CompareZScores(inputMatrix)
                        %CompareZScores() | Calculates all ZScores for matrix input
                        %	Calculates Zscores for each row of the input matrix
                        %	Inputs:
                        %		-inputMatrix{dataType} => 
                        %	Outputs:
                        %		-[zScoresVect]{dataType} => 
                        
                        sizeInputM = size(inputMatrix);
                        zScoresVect = zeros(sizeInputM(1),1);
                    end

                %Calc x from Z-score
                    function xVal = CalcXvalFromZscore(zScore,mu,stdDev)
                        xVal = zScore .* stdDev + mu;
                    end

                %Ztable Encode
					function [zTableM] = DisplayZTableValues()
					end

					function ProbDensityValue = FindProbDensityfromZtable(zScore)
					end

                %Probability Function
                    function output = CalcProb(x)
                        output = (1./sqrt(2.*pi())) .* exp(-1.*(x.^2./2));
                    end

                %Probability Density Function
                    function output = CalcProbDensity_zScore(zScore, varargin) %add precision via varargin 
                        syms x probFunc;
                        probFunc(x) = (1./sqrt(2.*pi)) .* exp(-1.*(x.^2 ./ 2));
                        %TRUE ProbDensityF VERSION:
                        %   output = integral([-infinity, zScore], probFunc(x))
                        %   - no way to calculate integrals on a computer, need to do it
                        %     numerically via trapz function
                        %   - no way to generate array to/from -infinity/+infinity,
                        %     therefore, need to exploit characteristics of the normal
                        %     distribution curve, like the fact that it sums to 1
                        %poop
                        
						% TODO: change from using trapz to the actual integral evaluated for speed, trapz takes too long

						% resolution declaration
							resolution = 0.0001;
							resolutionCounter = 4;
							if nargin > 1
								for i = 1:numel(varargin)
									if isfloat(varargin{i})
										resolution = varargin{i};
										resolutionCounter = 0;
										resolutionTemp = resolution;
										while resolutionTemp < 1
											resolutionCounter = resolutionCounter + 1;
											resolutionTemp = resolutionTemp .* 10.^resolutionCounter;
										end
										if resolution < 0.01
											resolutionCounter = resolutionCounter + 1;
										end
									end
								end
							end

                        output = 0;
                        
                        % Directional input
							% TODO: Add inputs to include >, >=, <, <=, and =
							direction = "left";
							for i=1:numel(varargin)
								if isstring(varargin{i})
									currentString = lower(varargin{i});
									if strcmp(lower(varargin{i}),"left") || strcmp(currentString,"<") || strcmp(currentString,"<=")
										direction = "left";
									elseif strcmp(lower(varargin{i}),"right") || strcmp(currentString,">") || strcmp(currentString,">=")
										direction = "right";
									elseif strcmp(lower(varargin{i}),"tocenter") || strcmp(currentString,"tomean")
										direction = "toCenter";
									end
								end
							end
						
						%Main math execution
						if strcmp(direction,"left") || strcmp(direction,"right")
							if (zScore == 0)
                            	output = 0.5;
                        	elseif (zScore > 0)
                            	Xs = 0:resolution:zScore; %make it so that if negative calculates from center out always then reverses to find values from left always
                            	Ys = probFunc(Xs);
								output = 0.5 + (trapz(Ys)./(10.^resolutionCounter)); %trapz area division has to be done according to how many samples on the x plane
							elseif (zScore < 0)
                            	Xs = zScore:resolution:0;
                            	Ys = probFunc(Xs);
                            	output = 0.5 - (trapz(Ys)./(10.^resolutionCounter)); %
                        	end

							if strcmp(direction,"right")
								output = 1 - output;
							end
						
						elseif strcmp(direction,"toCenter")
							if zScore > 0
								Xs = 0:resolution:zScore;
								output = (trapz(probFunc(Xs))./(10.^resolutionCounter)); %divide by size of samples
							elseif zScore < 0
								Xs = zScore:resolution:0;
								output = (trapz(probFunc(Xs))./(10.^resolutionCounter));
							end
						end

						%output section - determines if symbolic or not
							for i=1:numel(varargin)
								if isstring(varargin{i})
									if strcmp(varargin{i},"symbolic")
										return;
									end
								end
							end

							output = vpa(output);
                    end

					function output = CalcProbDensity_xVal(xVal, mu, sigma, varargin)
						zScore0 = StatsLib.CalcZScore(xVal, mu, sigma);
						
						if nargin > 3
							inputArguments = varargin{1:numel(varargin)};
							output = StatsLib.CalcProbDensity_zScore(zScore0, inputArguments);
						else
							output = StatsLib.CalcProbDensity_zScore(zScore0);
						end
					end

            % 5.2 - Solving Problems
                %ProbDensity between two z-Scores
                    function outputProbDensity = CalcProbDensityBetweenZScores(zScore1, zScore2)
                        container = [zScore1, zScore2];
						outputProbDensity = StatsLib.CalcProbDensity_zScore(max(container)) - StatsLib.CalcProbDensity_zScore(min(container));
                    end

                %ProbDensity between two x Values
                    function outputProbDensity = CalcProbDensityBetweenXVals(xVal1, xVal2, mu, sigma)
                        container = [StatsLib.CalcZScore(xVal1,mu,sigma), StatsLib.CalcZScore(xVal2,mu,sigma)];
						outputProbDensity = StatsLib.CalcProbDensityBetweenZScores(container(1), container(2));
                    end

				%ProbDensity for various xValues or zValues done in sequence
					function [outputProbDensitys] = CalcProbDensitys_FROM_InputTable(typeContainer, inputM)
						%FunctionName() short explanation here
						%	detailed description here
						%	Inputs:
						%		- typeContainer {vector<POINTERS>}          => declaration array outlining columnA types (zscore or xvals) and if xval, mu and sigma following first string
						%		- inputM        {vector<vector<POINTERS>>}  => table containing natural languange and operations to do in sequence to avoid typing the same command a bunch of times
						%			example: {Val a Direction ; val b Direction;...}
						%	Outputs:
						%		- output{dataType} => 
						
						if numel(typeContainer) == 1
							currentTypeDeclaration = lower(typeContainer(1));
							acceptableInputs = ["zscores" "zscore" "z scores" "z score"];
							if StatsLib.CheckUserStringInputValid(currentTypeDeclaration, acceptableInputs)
								sizeInputContainer = size(inputM);
								for i = 1:sizeInputContainer(1)
										
								end
							else
								error('typeContainer declaration incorrect');	
							end
						elseif numel(typeContainter) == 3
							currentTypeDeclaration = lower(typeContainer(1));
							if StatsLib.CheckUserStringInputValid(currentTypeDeclaration,acceptableInputs)
								
							else
								error('typeContainer type declaration incorrect');
							end
						else
							error('typeContainter input does not follow outlined input format')
						end
					end

            % 5.3 - Finding x values and z values from probability values
                %Percentile to zValue
                    function zScore = CalcZscoreFromPercentile(percentile, varargin)
                        zScore = Inf; %Assign fallback value that is impossible

                    end
            
                %Percentile to xValue
                    function xVal = CalcXvalFromPercentile(percentile, mu, sigma, varargin)
                        xVal = Inf; %set output value to something ridiculous for error checking
						%inverse ProbDensity functions to allow ProbDensity input to output X value
                    end
    end
end