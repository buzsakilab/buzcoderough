function AutoClustering(fbasename,elec,varargin)

% USAGE:
%     clu = AutoClustering(fbasename,elec)
% 
% AutoClustering automtically cleans the clu file defined by fbasename and
% electrode numbers. The program will look for the files
% fbasename.fet/res/spk.elec in the current folder.
% 
% INPUT:
%     fbasename: char array
%     elec: a vector of electrode numbers
%
% optional:
%     AutoClustering(fbasename,elec,dim)
%     where dim is the number of channels in electro group (if not
%     defined, will read the first line of the fet file
%     
% AutoClustering is meant to clean the output of KlustaKwik. The first
% thing it does is to separate electrical artifacts and MUA from putative
% isolated units. To do so, it sorts out units which have no clear
% refractory period (based on Hill, Mehta and Kleinfeld, J Neurosci.,
% 2012). Threshold can be set in the parameter section of this file
% ("Rogue spike threshold"). Then, it separates electrical
% artifats from MUA based on the assumption that electrical artifacts are
% highly correlated on the different channels: the average waveform of at
% least one channel has to be different from the across-channel average
% waveform by a certrain amount of total variance (can be set in the
% parameter section, "Deviation from average spike threshold")
%
%
% Once the program has determined which of the clusters are putative
% isolated units, it tries to merge them based on waveform similarity
% (mahalanobis distance) and quality of the refractory period in the new
% merged cluster (or "Inter Common Spike Interval" from MS Fee et al. J 
% Neurosci. Meth., 1996)
%
% Adrien Peyrache, 2012



% Parameters
% Number recording sites
if ~isempty(varargin)
    dim = varargin{1};
    dim = dim(:);
    if any(double(int16(dim))~=dim)
        error('Number of dimensions must be an integer')
    end
      
    if size(dim,1) ~= length(elec) && length(dim) ~=1
        error('Number of dimensions must be a vector of the same length as electrode vecotr or a single value')
    end
    if length(dim) == 1
        dim = dim*ones(length(elec),1);
    end
else
    dim = zeros(length(elec),1);
end

% Load spike waveforms? Used for detection of electrical artifacts
loadspk = 0;
% Refractory period in msec
tR = 1.5;
% Censored period in msec (specific to the spike detection process)
tC = 0.85;
% Rogue spike threshold (for MUA); value between 0 an 1
%rogThres = 0.25;
rogThres = 0.25;
% Relative deviation (from 0 to 1) from average spike threshold (for electrical artifacts)
%devThres = 0.25;
devThres = 1000;

% Do Merging?
doMerge = 0;

% overwrite clu file
rewriteclu= 1;
% Write a log file?
WriteLogFile = 1;

if WriteLogFile
    tic;
    log = [];
end

elec = elec(:)';
if length(elec)>1
    for eIx=1:length(elec)
        AutoClustering(fbasename,elec(eIx),dim(eIx));
    end
    
else

    % Load fet, clu and res
    fprintf('Sorting electrode %i of %s\n',elec,fbasename)
    fet = LoadFeatures(fbasename,elec,dim);

    clu = load([fbasename '.clu.' num2str(elec)]);
    clu = clu(2:end);
    res = load([fbasename '.res.' num2str(elec)])/20;

    % load waveforms
    if loadspk
        %wav = LoadSpikeWaveforms([fbasename '.spk.' num2str(elec)],dim,32);
        
    end

    % number of spikes
    n = size(fet,1);

    %Power of each waveform - UNUSED so far
    %pw = sqrt(sum(fet.^2,2));
    % Percentile of the power
    % p = percentile(pw,0.99);

    %  devFromMean quantifies how much the spike are different on
    %  the different channels (what is the maximal distance from the averaged
    %  spike)
    devFromMean = [];

    fractRogue = [];
    meanPw = [];

    for ii=1:length(unique(clu))
        rg = res(clu==ii);
        
        if loadspk
            m = squeeze(mean(wav(:,:,clu==ii),3));
            y = pdist(m,'euclidean')/max(sqrt(sum(m.^2,2)));
            devFromMean = [devFromMean;max(y)];
        end

        l = FractionRogueSpk(rg,tR,tC);    
        fractRogue = [fractRogue;l];

        %meanPw = [meanPw;sum(pw(clu==ii)>p)/sum(clu==ii)];   

    end

    % Definition of cluster 0 (noiseIx) and cluster 1 (muaIx) 
    % Outliers of total spike power (putative electrical artifacts) not implemented yet

    if loadspk
        noiseIx = find(fractRogngue>rogThres & devFromMean<devThres);
        muaIx = find(fractRogue>rogThres & devFromMean>devThres);
    else
        muaIx = find(fractRogue>rogThres);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display intermediary results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if 0
        for ii=1:length(unique(clu))
            rg = res(clu==ii);
            [ac,b] = CrossCorr(rg*10,rg*10,1,60);ac(b==0)=0;
            figure(1),clf
                bar(b,ac,1)
                disp([fractRogue(ii) devFromMean(ii)])
                pause
        end
        for ii=1:length(muaIx)
            rg = res(clu==muaIx(ii));
            [ac,b] = CrossCorr(rg*10,rg*10,1,60);ac(b==0)=0;
            figure(1),clf
                bar(b,ac,1)
                disp([fractRogue(muaIx(ii)) devFromMean(muaIx(ii))])
                pause
        end
    end

    % Here we compute # of spike per cell. Some code for the errormatrix fails
    % when the cluster is defined by only a few samples. We'll put a
    % threshopld a bit later on the total # of spikes.
    h = hist(clu,unique(clu));
    h = h(:);

    % Update the clu indices
    if loadspk
        for ii=1:length(noiseIx)
            clu = updateclu(clu,noiseIx(ii),1000+ii);
            log = [log sprintf('%d -> %d; Looks like electrical artifact\n',noiseIx(ii),1000+ii)];
        end
    end
    for ii=1:length(muaIx)
        clu = updateclu(clu,muaIx(ii),2000+ii);
        log = [log sprintf('%d -> %d; Looks like MUA\n',muaIx(ii),2000+ii)];
    end

    % Here we select only clusters that correspond to putative units and that
    % have at least 20 spikes (otherwise errormatrix calculation fails)
    goodCluIx = ismember(clu,find(fractRogue<rogThres & h>20));

    % if there is anything to compare...
    mergehistory = [];
 
    if length(unique(clu(goodCluIx)))>1
        if max(clu(goodCluIx)) ~= max(clu(clu<1000))
                ix = clu == max(clu(goodCluIx));
                clu(ix) = max(clu(clu<1000))+1;
        end
        newclu = clu(goodCluIx);
        
        if doMerge
            % merge similar clusters which are neither noise nor MUA
            try
                [newclu mergehistory] = mergeclu_slow(mewclu,res(goodCluIx),fet(goodCluIx,:),tR,tC,doMerge);
                % reassign cluster indices
               
            catch
                warning(['merging failed: ' lasterr])
            end
        else
            
            em = errormatrix(fet(goodCluIx,:),newclu);
            ems = max(em,em');   
            y = squareform((1-ems)-diag(diag(1-ems)),'tovector');
            Z = linkage(y);
            [T,perm] = dendrogram_struct(Z,0);
            
            cluIx = unique(newclu);
            for ii=1:length(cluIx)
                newclu = updateclu(newclu,cluIx(perm(ii)),max(cluIx)+ii+1);
            end
        end
        
        clu(goodCluIx) = newclu;        

    end
    
    cluIx = unique(clu(clu<1000));
    clunew = clu;
    for ii=1:length(cluIx)
        clunew = updateclu(clunew,cluIx(ii),ii+1);
    end
    clu = clunew;
    
    % Log
    for ii=1:size(mergehistory,1)
        log = [log sprintf('%d + %d -> %d\n',mergehistory(ii,1),mergehistory(ii,2),mergehistory(ii,3))];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display final results (requires function CrossCorr, not in the
    % toolbox)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if 0
    cluIx = unique(newclu(newclu>1));
    for x=1:length(cluIx)
        rgx = res(newclu==cluIx(x));
        [ax,b] =  CrossCorr(rgx*10,rgx*10,1,60);ax(b==0)=0;
        for y=x+1:length(cluIx)
            rgy = res(newclu==cluIx(y));
            [ay,b] =  CrossCorr(rgy*10,rgy*10,1,60);ay(b==0)=0;
            [fxy fcorr] = crossrefract(rgx,rgy);
            [h,b] = CrossCorr(rgx*10,rgy*10,1,60);ac(b==0)=0;
            figure(1),clf
            subplot(1,3,1)
                bar(b,ax,1)
             subplot(1,3,2)
                bar(b,h,1)
                title([fxy fcorr])
             subplot(1,3,3)
                bar(b,ay,1)
            pause
        end
    end
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Write new clu file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rewriteclu

        fid = fopen([fbasename '.clu.' num2str(elec)],'w');
        fprintf(fid,'%i\n',[length(unique(clu));clu]);
        fclose(fid);
    end

    if WriteLogFile
        % Create (of overwrite) a log file
        fid = fopen([fbasename '.alg.' num2str(elec)],'w');
        time_tot = toc;
        fprintf(fid,[log 'That took %f seconds\n'],time_tot);
    end
end
