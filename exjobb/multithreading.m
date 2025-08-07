function multithreading()
% Use Parallel pools to speed up the code.
    myCluster = parcluster('local');
    NumWorkers = myCluster.NumWorkers;
    p = gcp("nocreate");
    if isempty(p) % If no pool, create new one.
        parpool(NumWorkers);
    else % If there is a pool, but its small, delete it and create new one.
        if p.NumWorkers < NumWorkers
            delete(gcp('nocreate'));
            parpool(NumWorkers);
        end
    end
end