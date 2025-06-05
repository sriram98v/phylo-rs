#!/bin/bash

ITERATIONS=1000

# timeout limit in seconds
TIME_OUT=3

TREE_SIZE_START=200
TREE_SIZE_STEP=200
TREE_SIZE_END=10000

OUTPUT_FILE=runtimes.csv

data1='./data/t1/'
data2='./data/t2/'

mkdir -p ${data1}
mkdir -p ${data2}

echo "Creating trees..."
for i in $(seq $TREE_SIZE_START $TREE_SIZE_STEP $TREE_SIZE_END); do
  gotree generate yuletree --seed 10 -l ${i} -o ${data1}/${i}.nwk
  gotree generate yuletree --seed 12 -l ${i} -o ${data2}/${i}.nwk
done


function run_tests() {

    SCRIPT=${1}

    SCRIPT_NAME=$(basename "${1}")

    ALGO_NAME=${SCRIPT_NAME%.*}

    DATA_DIR=${2}
    OUTPUT=${3}

    echo -en ${ALGO_NAME} >> ${OUTPUT}

    printf "%-40s\n" ${SCRIPT_NAME}
    for ntaxa in $(seq $TREE_SIZE_START $TREE_SIZE_STEP $TREE_SIZE_END); do
        local min=0 max=0 sum=0 avg
        
        for i in $(seq 1 ${ITERATIONS}); do


            script_out=`timeout ${TIME_OUT}s /usr/bin/time -v ${SCRIPT} ${DATA_DIR}/t1/${ntaxa}.nwk ${DATA_DIR}/t2/${ntaxa}.nwk 2>&1`

            success=$?

            # echo $success
            # Break if the test failed. We do not need to repeat it then.
            if [[ ${success} != 0 ]] && [[ ${success} != 124 ]]; then
                echo ${success}
                echo ${script_out}
                break
            fi

            duration=`echo ${script_out} | sed "s/.*Internal time: \([0-9.e-]*\).*/\1/g" | sed "s/e/\\*10\\^/" | sed "s/+//"`
            duration=`echo "scale=7; ${duration}" | bc`

            if [[ ${success} == 124 ]]; then
                duration=${TIME_OUT}
                # min=${duration}
                # max=${duration}
                # avg=${duration}
                # echo -en ",${avg}" >> ${OUTPUT}
                # break
            fi

            if [[ ${max} == 0 ]]; then
                min=${duration}
                max=${duration}
            else
                if [[ `echo "${duration} > ${max}" | bc` == 1 ]]; then
                    max=${duration}
                fi
                if [[ `echo "${duration} < ${min}" | bc` == 1 ]]; then
                    min=${duration}
                fi
            fi
            sum=`echo "${sum} + ${duration}" | bc`

        done

        if [[ ${success} == 0 ]] || [[ ${success} == 124 ]]; then
            if [[ ${success} == 0 ]]; then
                avg=`echo "scale=7;${sum}/${ITERATIONS}" | bc`
            fi

            printf "% 15s " ${ntaxa}
            printf "% 5.7fs " ${min}
            printf "% 5.7fs " ${max}
            printf "% 5.7fs\n" ${avg}


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

    echo "Using trees in $DATA_DIR"
    
    echo "Writing runtimes to $OUPUT_FILE"


    for algo in ${ALGOS_DIR}/*;do
        if [ -x "$algo" ]; then
            run_tests $algo ${DATA_DIR} ${OUPUT_FILE}
        fi
    done

}

echo "Deleting ${OUTPUT_FILE}"
rm ${OUTPUT_FILE} 2> /dev/null

echo -en "algorithms" >> ${OUTPUT_FILE}
for ntaxa in $(seq $TREE_SIZE_START $TREE_SIZE_STEP $TREE_SIZE_END); do
    echo -en ",$ntaxa" >> ${OUTPUT_FILE}
done
echo -en "\n" >> ${OUTPUT_FILE}


run_all_tests phylotree/time data runtimes.csv
run_all_tests compacttree/time data runtimes.csv
run_all_tests genesis/time data runtimes.csv
run_all_tests treeswift/time data runtimes.csv
run_all_tests dendropy/time data runtimes.csv
run_all_tests phylo-rs/time data runtimes.csv
run_all_tests ape/time data runtimes.csv