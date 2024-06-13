#! /bin/bash
docker run --rm  --volume "$PWD/v2024.1:/data"  --user "${UID}:${GROUPS[0]}"  --env JOURNAL=joss  openjournals/paperdraft
