function A = mazeFineSpectrum(A)

A = getResource(A,'PfcTrace');
A = getResource(A,'Sleep1Epoch');
A = getResource(A,'Sleep2Epoch');
sleepEpoch{1} = sleep1Epoch{1};
sleepEpoch{2} = sleep2Epoch{1};

A = registerResource(A, 'SleepLargeSpecgram', 'tsdArray', {[],[]}, ...
    'sleepLargeSpecgram', ...
    'sleepclarge spectrum','mfile');

A = registerResource(A, 'SleepLargeSpecgramFreq', 'numeric',  {[],[]}, ...
    'sleepLargeSpecgramFreq', ...
    ['frequencies [0-200 Hz] for the sleep specgrams'],'mfile');



[p,dset,e] = fileparts(current_dir(A));

eegfname = [current_dir(A) filesep dset 'eeg' num2str(pfcTrace) '.mat'];
if exist([eegfname '.gz'])
    display(['unzipping file ' eegfname]);
    eval(['!gunzip ' eegfname '.gz']);
end
load(eegfname)
display(['zipping file ' eegfname]);
%  eval(['!gzip ' eegfname]);

eval(['eegPfc = EEG' num2str(pfcTrace) ';']);
eval(['clear EEG' num2str(pfcTrace)]);

params.Fs = 2083;
params.fpass = [20 100];
params.err = [2, 0.95];
params.trialave = 0;
movingwin = [0.2, 0.1];

for s=1:2
	
	eeg = Restrict(eegPfc, sleepEpoch{s});
	
	%  keyboard
	
	Fs = floor(1 / median(diff(Range(eeg, 's'))));
		
	dEeg = Data(eeg);
	[S,t,f,Serr]=mtspecgramc(dEeg,movingwin,params);
		
	%all this fuss is necessary to accommodate for EEG recordigns with possible
	%gaps in them
		times = Range(eeg);
		t1 = 0:(1/Fs):((1/Fs)*(length(times)-1));
		[t2, ix] = Restrict(ts(t1), t-movingwin(1)/2);
		times = times(ix);
		
		sleepLargeSpecgram{s} = tsd(times, S);
	
	display(['sleep' num2str(s) 'ok']);

end

%  figure(1),clf
%  imagesc(Range(mazeFineSpecgramPfc{1}),f,log10(Data(mazeFineSpecgramPfc{1})+eps)');
%  title('Sleep1')
%  axis xy
%  
%  figure(2),clf
%  imagesc(Range(mazeFineSpecgramPfc{2}),f,log10(Data(mazeFineSpecgramPfc{2})+eps)');
%  title('Sleep1')
%  axis xy


sleepLargeSpecgramFreq = f;
sleepLargeSpecgram = tsdArray(sleepLargeSpecgram);

A = saveAllResources(A);

