function MCparforCollProb(numWorkers)
matlabpool('local', numWorkers)

% user defined parallel operations must occur inside an spmd block
spmd
    fprintf(['hello! I am ' num2str(labindex) ' of ' num2str(numlabs) '\n']);
end

ProbcumMC1000=cell(1000,1);
parfor i=1:1:1000
ProbcumMC1000{i}=MC_cummCollProb(i);
end

save('MCcummProbColl_1000','ProbcumMC1000')








% all done, close pool of workers
matlabpool close
exit;