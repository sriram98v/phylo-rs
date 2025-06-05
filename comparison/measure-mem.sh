#!/bin/bash

ITERATIONS=10

OUTPUT_FILE=memory-utilization.csv

data='./mem-trees/'

mkdir -p ${data}

echo "Creating trees..."
for i in 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000 ; do
  gotree generate yuletree --seed 10 -l ${i} -o ${data}/${i}.nwk
done


function run_tests() {

    SCRIPT=${1}

    SCRIPT_NAME=$(basename "${1}")

    ALGO_NAME=${SCRIPT_NAME%.*}

    DATA_DIR=${2}
    OUTPUT=${3}

    echo -en ${ALGO_NAME} >> ${OUTPUT}

    printf "%-40s\n" ${SCRIPT_NAME}
	for ntaxa in 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000 ; do
		local min=0 max=0 sum=0 avg

        for i in $(seq 1 ${ITERATIONS}); do


            script_out=`/usr/bin/time -v ${SCRIPT} ${DATA_DIR}/${ntaxa}.nwk 2>&1`

            success=$?


            # echo ${success}
            # echo ${script_out}

            # echo $success
            # Break if the test failed. We do not need to repeat it then.
            if [[ ${success} != 0 ]]; then
                echo ${success}
                echo ${script_out}
                break
            fi

            mem=`echo ${script_out} | sed "s/.*Maximum resident set size .kbytes.: \([0-9]*\).*/\1/g"`
			mem=`echo "scale=3;${mem}/1024" | bc`


            if [[ ${max} == 0 ]]; then
                min=${mem}
                max=${mem}
            else
                if [[ `echo "${mem} > ${max}" | bc` == 1 ]]; then
                    max=${mem}
                fi
                if [[ `echo "${mem} < ${min}" | bc` == 1 ]]; then
                    min=${mem}
                fi
            fi
            sum=`echo "${sum} + ${mem}" | bc`

        done


        if [[ ${success} == 0 ]]; then
            if [[ ${success} == 0 ]]; then
                avg=`echo "scale=3;${sum}/${ITERATIONS}" | bc`
            fi

            printf "% 15s " ${ntaxa}
            printf "% 5.3f Mb " ${min}
            printf "% 5.3f Mb " ${max}
            printf "% 5.3f Mb\n" ${avg}


            # Print to tab files for easier post-processing
            echo -en ",${avg}" >> ${OUTPUT}
        else
            echo "Fail!"
        fi
    done
    echo -en "\n" >> ${OUTPUT}

    return ${success}
}

function run_all_tests() {
    ALGOS_DIR=${1}

    DATA_DIR=${2}
    OUPUT_FILE=${3}

    echo $ALGOS_DIR
    
    echo "Using trees in $DATA_DIR"
    
    echo "Writing runtimes to $OUPUT_FILE"


    for algo in ${ALGOS_DIR}/*;do
        echo $algo
        if [ -x "$algo" ]; then
            run_tests $algo ${DATA_DIR} ${OUPUT_FILE}
        fi
    done

}

echo "Deleting ${OUTPUT_FILE}"
rm ${OUTPUT_FILE} 2> /dev/null

echo -en "algorithms" >> ${OUTPUT_FILE}
for ntaxa in 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000 ; do
    echo -en ",$ntaxa" >> ${OUTPUT_FILE}
done
echo -en "\n" >> ${OUTPUT_FILE}

run_all_tests phylotree/mem mem-trees ${OUTPUT_FILE}
run_all_tests compacttree/mem mem-trees ${OUTPUT_FILE}
run_all_tests genesis/mem mem-trees ${OUTPUT_FILE}
run_all_tests treeswift/mem mem-trees ${OUTPUT_FILE}
run_all_tests dendropy/mem mem-trees ${OUTPUT_FILE}
run_all_tests phylo-rs/mem mem-trees ${OUTPUT_FILE}
run_all_tests ape/mem mem-trees ${OUTPUT_FILE}
