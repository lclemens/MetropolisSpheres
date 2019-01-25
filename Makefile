all:
	gcc -O3 driveMetropolis.c -o metropolis.out -lm


run:
	gcc -O3 driveMetropolis.c -o metropolis.out -lm

	./metropolis.out parameters.txt output 6

batch:

	gcc -O3 driveMetropolis.c -o metropolis.out -lm

	chmod a+x runBatch.sh

	./runBatch.sh

hpc:

	gcc -O3 driveMetropolis.c -o metropolis.out -lm
	
	chmod a+x makePubs.pl

	./makePubs.pl

	chmod a+x submit.pl

	./submit.pl

copyrepo:

	git -C ~/Documents/PolymerGit log -1 --pretty=format:%H > "CommitUsedHash.txt"

	cp -r ~/Documents/PolymerGit/src/PolymerCode/ .

parallel:

	gcc -O3 driveMetropolis.c -o metropolis.out -lm

	chmod a+x batchMetropolis.sh

	./batchMetropolis.sh 7
