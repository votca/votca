docker run --rm  --volume $PWD/v2024:/data  --user $(id -u):$(id -g)  --env JOURNAL=joss  openjournals/paperdraft
