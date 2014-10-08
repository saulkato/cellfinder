%Count cells from static Z stack demo
%
%use Laplacian of Gaussian to find blobs with some extra tricks
%use BlobThread heuristic to bind blobs
%
%Saul Kato

%% beads

stackFile1='beads_177.32 dwell_lastone.tif';   %this is an 8-bit stack, watch out

findBlobOptions.thresholdMargin=30;  %low because it is 8-bit
findBlobOptions.minSpacing=2;
findBlobOptions.shapeFilterWidth=2;
findBlobOptions.medianFilterWidth=1;

bindBlobOptions.maxBindingDistance=2;
bindBlobOptions.minThreadLength=2;

[cells, blobs]=FindCellsStaticStack(stackFile1,findBlobOptions,bindBlobOptions);

disp(['Found ' num2str(cells.n) ' cells in ' stackFile1 '.']);



%% mouse cells

stackFile2='C1-ZSeries-08072014-1700-001_Cycle00001_CurrentSettings_RedChannel.tif';

findBlobOptions.thresholdMargin=500;  %low because it is 8-bit
findBlobOptions.minSpacing=2;
findBlobOptions.shapeFilterWidth=8;
findBlobOptions.medianFilterWidth=1;

bindBlobOptions.maxBindingDistance=4;
blobBondOptions.minNumZplanes=4;


[cells, blobs]=FindCellsStaticStack(stackFile2,findBlobOptions,bindBlobOptions);

disp(['Found ' num2str(cells.n) ' cells in ' stackFile2 '.']);
