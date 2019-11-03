%% PROJECT CONFIG

PROJECT_DIR = pwd ; % put base path to project her
cd(PROJECT_DIR)

%% add to the path

projPathDirs = { 
    'src' 
    'data'
    'bin'
} ;

for idx=1:length(projPathDirs)
    addpath(genpath(strcat(PROJECT_DIR,'/',projPathDirs{idx})))
end

clear projPathDirs

%% setup global vars

OUTSTR = 'run1' ;

%% make output directory vars

OUTDIR = strcat(PROJECT_DIR , '/data/') ; 
OUTDIR_INTERM = strcat(OUTDIR, '/interim/' ) ;
OUTDIR_PROC = strcat(OUTDIR, '/processed/' ) ;

%%

PARC_STR = 'schaefer200-yeo17' ;
NUM_NODES = 200 ; 

COMM_PROB_THR = 0.05 ;

