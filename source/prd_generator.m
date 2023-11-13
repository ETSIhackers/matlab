function prd_generator(fullFileName)
% PRD_GENERATOR PET raw data generator
% !! THIS DOES NOT WRITE VALID H5 FILES
%
% See https://github.com/ETSInitiative/PRDdefinition
%
% Note detector ids, energy bin indices and tof indices start at 0
% 
% Example:
%  ffn = tempname ;
%  prd_generator(ffn)
%  prd_analysis(ffn)
%
% See also prd_analysis


% detector id Starts at 0 in this MATLAB code
% Format of hdf5 file unclear from GitHub code - revese engineer from h5
% example

arguments
    fullFileName 
end

% // these are constants for now (from PRDdefinition GitHub 25 Oct 2023)
NUMBER_OF_ENERGY_BINS = 3;
NUMBER_OF_TOF_BINS = 300;
RADIUS = 400. ;
NUMBER_OF_TIME_BLOCKS = 6;
NUMBER_OF_EVENTS = 1000 ;   % this one is missing from the Python ('cos its different for each time block?)
COUNT_RATE = 500 ;
tofResolution = 9.4 ;
energyResolutionAt511 = 0.11 ;
listmodeTimeBlockDuration = 1 ;

lowestEnergyBinEdge  = 430 ; %keV
highestEnergyBinEdge = 650 ;

number_of_detectors = 10 ; 
nangles = number_of_detectors ;  % Assume equal angular spacing


% What does x,y,z mean in terms of patient or gantry coordinates?


for iangle = 1: nangles
    angle_rad = (iangle-1)*2*pi / nangles ; % detector angular position in radians

    x(iangle) = RADIUS * sin(angle_rad) ;
    y(iangle) = RADIUS * cos(angle_rad) ;
    z(iangle) = 0 ;
    id(iangle) = iangle - 1 ; % detector ids start at 0.
end
detectors.x = x; detectors.y = y; detectors.z = z ; detectors.id = id ;


% TOF Bin Info (mm
% !! Bin edges not currently used. Appears to be from -radius to radius in NUMBER_OF_TOF_BINS+1 edges, 
% but might be confusing as to how the array index maps to edge? Better to
% have an array which is [nbin 2]?
% Needs to specify a direction, or how it corresponds to angle?
tofBinEdges = linspace(-RADIUS, RADIUS, NUMBER_OF_TOF_BINS+1)' ;

% Energy Bins. Assumed equal spaced, incremental order and no gaps (keV?)
% ! Not actually written out

energyBinEdges = linspace(lowestEnergyBinEdge, highestEnergyBinEdge, NUMBER_OF_ENERGY_BINS+1) ;
energyBinEdges = energyBinEdges(:) ;


% Coincendence events. Scanner cannot distinguish source of coincidence
% (true, random, scatter etc) so although these may be in the simulation,
% all that is written to file is coincidence events.

% There are NUMBER_OF_TIME_BLOCKS time blocks, and within each the number
% of events follows a Poisson distribution. Here we pre-determine the
% number of events per time block. OK for simulation, but in reality we
% cannot know ahead of time. Issue is in allocating space and writing to h5
% file when total counts not known until the end of the scan.

time_block_ids    = zeros([NUMBER_OF_TIME_BLOCKS 1]) ;
time_block_counts = zeros([NUMBER_OF_TIME_BLOCKS 1]) ;

for itblock = 1: NUMBER_OF_TIME_BLOCKS
    time_block_ids(itblock) = itblock - 1 ;
    time_block_counts(itblock) = poissrnd(COUNT_RATE) ;
end

num_events_total = sum(time_block_counts(:)) ;

eventDatatype = 'uint32' ; 
eventSource   = 'positron-annihilation' ;
   

dtb.id = uint32(time_block_ids) ;
promptEvents = cell([]) ;

% % for itime_block = 1: NUMBER_OF_TIME_BLOCKS
% %     promptEvents{itime_block} = get_events(time_block_counts(itime_block), eventSource , ...
% %           number_of_detectors, NUMBER_OF_ENERGY_BINS, NUMBER_OF_TOF_BINS, eventDatatype) ;
% % end
% % dtb.promptEvents = promptEvents ;


% H5 writing
% Header scanner and exam
dh.scanner.modelName.value = {''} ;
dh.scanner.modelName.has_value = uint32(0) ;

dh.scanner.detectors = {detectors} ;

dh.scanner.tofBinEdges = {single(tofBinEdges)} ;
dh.scanner.tofResolution = single(tofResolution) ;
dh.scanner.energyBinEdges = {single(energyBinEdges)} ;
dh.scanner.energyResolutionAt511 = single(energyResolutionAt511) ;
dh.scanner.listmodeTimeBlockDuration = single(listmodeTimeBlockDuration) ;

dh.exam.value.institution.name = {'Diamond Light Source'} ;
dh.exam.value.institution.address = {'Harwell Science and Innovation Campus, Didcot, Oxfordshire, OX11 0DE, UK'} ;

dh.exam.value.subject.id = {'123456'} ;


% Write header to h5 file

h5create(fullFileName, "/PrdExperiment/header/scanner/modelName/value", [1 1], 'Datatype', 'string' )
 h5write(fullFileName, "/PrdExperiment/header/scanner/modelName/value", dh.scanner.modelName.value)

% % Confusing what PRD are doing here with nested datasets?
% % Confusing about capitals in dataset names vs variables
% 
% h5create(fullFileName,"/Header/ExamInformation/Subject/id", [1 1], 'Datatype', 'string')
%  h5write(fullFileName,"/Header/ExamInformation/Subject/id", subject.id)
% 
% h5create(fullFileName,"/Header/ExamInformation/Institution/name", [1 1], 'Datatype', 'string') ;
%  h5write(fullFileName,"/Header/ExamInformation/Institution/name", institution.name) ;
% 
% h5create(fullFileName,"/Header/ExamInformation/Institution/address", [1 1], 'Datatype', 'string') ;
%  h5write(fullFileName,"/Header/ExamInformation/Institution/address", institution.address) ;
% 
% % For detectors, PRD might be using HDF5 composite (compound) variables
% % Pain to do that here without low-level hdf5 functions
% 
% h5create(fullFileName,"/Header/ScannerInformation/Detector", [nangles 4]) ;
%  h5write(fullFileName,"/Header/ScannerInformation/Detector",...
%     [detectors.x ; detectors.y ; detectors.z ; detectors.id]')  ;
% 
% % Write event data 
% h5create(fullFileName,"/TimeBlock/TimeBlockId", [1 NUMBER_OF_TIME_BLOCKS]) ;
%  h5write(fullFileName,"/TimeBlock/TimeBlockId", time_block_ids) ;
% 
%  h5create(fullFileName,"/TimeBlock/TimeBlockCount", [1 NUMBER_OF_TIME_BLOCKS]) ;
%   h5write(fullFileName,"/TimeBlock/TimeBlockCount", time_block_counts) ;
% 
% h5create(fullFileName,"/TimeBlock/CoincidenceEvent", [num_events_total 5], 'Datatype', eventDatatype) ;
%  h5write(fullFileName,"/TimeBlock/CoincidenceEvent", ...
%      [coincidenceEvents.detector_1_id ; ...
%       coincidenceEvents.detector_2_id ; ...
%       coincidenceEvents.energy_1_idx ; ...
%       coincidenceEvents.energy_2_idx ; ...
%       coincidenceEvents.tof_idx ]') ;

h5disp(fullFileName) 

end

% - - - - - 

function coincidenceEvents = get_events(nevents, eventSource, number_of_detectors, NUMBER_OF_ENERGY_BINS, NUMBER_OF_TOF_BINS, eventDatatype)
% coincidenceEvents 
% coincidenceEvents = get_events(nevents, eventSource, number_of_detectors, NUMBER_OF_ENERGY_BINS, NUMBER_OF_TOF_BINS, eventdatatype)
%

% I think I may have mis-understood and we just have positron annihilation 
switch eventSource
    case {'positron-annihilation'}
        % Coincidence Events
        coincidenceEvents = struct([]) ;  % May be better way to handle these. Need to understand hdf5 layout first

        for ievent = 1: nevents
            % Following are 0-based
            % randn generates numbers evenly spaced 0 to 1. 
            % We are going to have to use 0-based indexing (for compatibility), hence use of floor: 
            coincidenceEvents.detector1id(ievent) = cast( floor( rand * number_of_detectors )  , eventDatatype) ;
            coincidenceEvents.detector2id(ievent) = cast( floor( rand * number_of_detectors )  , eventDatatype);
            coincidenceEvents.energy1idx(ievent)  = cast( floor( rand * NUMBER_OF_ENERGY_BINS) , eventDatatype) ;
            coincidenceEvents.energy2idx(ievent)  = cast( floor( rand * NUMBER_OF_ENERGY_BINS) , eventDatatype) ;
            coincidenceEvents.tofidx(ievent)      = cast( floor( rand * NUMBER_OF_TOF_BINS )   , eventDatatype);
        end
    otherwise
        error(['Unknown eventSource: ',eventSource])
end


end

