docker run --rm  --volume $PWD/v2022.1:/data  --user $(id -u):$(id -g)  --env JOURNAL=joss  openjournals/paperdraft
