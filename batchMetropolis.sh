# Parallelization Script
# July 13, 2016

    # set number of runs submitted by checking the running processes for lines with the program name, set NRequested to number of lines (may include one more than actual number of runs, since it counts grep -c metropolis in its tally)

NRequested=`ps | grep -c metropolis`

# while number of iterations ran is less than or equal to total number of iterations desired, loop through runs

ITERATIONS=1

for((SRADIUS=1;SRADIUS<11;SRADIUS++))
do


    # loop to periodically check how many runs are submitted

    while (( $NRequested >= $1 ))   # while number requested is greater than number of processors we want to use (user input in command line)

        do
            sleep 1     # wait 1 sec before checking again

            NRequested=`ps | grep -c metropolis` # check again

    done

        echo "Done sleeping."

    # run program with specified parameters

    ./metropolis.out parameters.txt sphereLocations.$SRADIUS $SRADIUS &

    # update number of running programs
    NRequested=`ps | grep -c metropolis`

    echo "Done calling metropolis."
done

# wait for all background processes to finish before concatenating files
wait

echo "Done waiting for processes to finish."


