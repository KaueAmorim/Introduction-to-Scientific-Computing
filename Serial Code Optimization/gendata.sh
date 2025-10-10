#!/bin/bash

PROG=matmult
CPU=3

DATA_DIR="Dados/"

mkdir -p ${DATA_DIR}

echo "performance" > /sys/devices/system/cpu/cpufreq/policy${CPU}/scaling_governor

make purge matmult

METRICA="FLOPS_DP L3CACHE ENERGY"
TAMANHOS="64 100 128 1024 2000"
TEMPOS="${DATA_DIR}/Tempos.csv"

for m in ${METRICA}
do
    LIKWID_CSV="${DATA_DIR}/${m}.csv"
    rm -f ${TEMPOS}
    rm -f ${LIKWID_CSV}

    for n in $TAMANHOS
    do
        LIKWID_OUT="${DATA_DIR}/${m}_${n}.txt"
        
        echo "--->>  $m: ./${PROG} $n" >/dev/tty
        # Executa com LIKWID
        likwid-perfctr -O -C ${CPU} -g ${m} -o ${LIKWID_OUT} -m ./${PROG} ${n} >> ${TEMPOS}
    done

    # Define o nome da métrica a ser procurada com base no grupo
    case "$m" in
        "FLOPS_DP")
            VARIAVEL="FP_ARITH_INST_RETIRED_SCALAR_DOUBLE"
            POS=3
        ;;
        "L3CACHE")
            VARIAVEL="L3 miss ratio"
            POS=2
        ;;
        "ENERGY")
            VARIAVEL="Energy \[J\]"
            POS=2
        ;;
    esac

    LIKWID_CSV="${DATA_DIR}/${m}.csv"

    for n in $TAMANHOS
    do
        LIKWID_OUT="${DATA_DIR}/${m}_${n}.txt"
        # Escreve linha no CSV
        echo "${n},$(grep -E "${VARIAVEL}" ${LIKWID_OUT} | cut -d ',' -f${POS} | xargs | tr ' ' ',')" >> ${LIKWID_CSV}
    done
done

echo "powersave" > /sys/devices/system/cpu/cpufreq/policy${CPU}/scaling_governor