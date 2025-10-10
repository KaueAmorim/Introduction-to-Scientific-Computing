#!/bin/bash

PROG=matmult
CPU=3

DATA_DIR="Dados/"

mkdir -p ${DATA_DIR}

echo "performance" > /sys/devices/system/cpu/cpufreq/policy${CPU}/scaling_governor

make purge matmult

METRICA="L3CACHE ENERGY FLOPS_DP"
TAMANHOS="64 100 128 1024 2000 2048 3000 4096 6000 7000 10000 50000 60000 70000 100000"
TEMPOS="${DATA_DIR}/Tempos.csv"

# Documenta a arquitetura do processador
echo "Documentando arquitetura do processador..." >/dev/tty
likwid-topology -g -c > ${DATA_DIR}/arquitetura_processador.txt

rm -f ${TEMPOS}
echo "N,t_matVet,t_matVet_otim,t_matMat,t_matmat_otim" > ${TEMPOS}

for m in ${METRICA}
do
    LIKWID_CSV="${DATA_DIR}/${m}.csv"
    rm -f ${LIKWID_CSV}
    
    # Define cabeçalho baseado na métrica
    case ${m} in
        "L3CACHE")
            echo "N,L3_miss_ratio_matVet,L3_miss_ratio_matVet_otim,L3_miss_ratio_matMat,L3_miss_ratio_matmat_otim" > ${LIKWID_CSV}
            METRIC_PATTERN="L3 miss ratio"
            ;;
        "ENERGY")
            echo "N,Energy_J_matVet,Energy_J_matVet_otim,Energy_J_matMat,Energy_J_matmat_otim" > ${LIKWID_CSV}
            METRIC_PATTERN="Energy \[J\]"
            ;;
        "FLOPS_DP")
            echo "N,MFLOPS_matVet,MFLOPS_matVet_otim,MFLOPS_matMat,MFLOPS_matmat_otim" > ${LIKWID_CSV}
            METRIC_PATTERN="MFLOP/s"
            ;;
    esac

    for n in $TAMANHOS
    do
        LIKWID_OUT="${DATA_DIR}/${m}_${n}.txt"
        
        echo "--->>  $m: ./${PROG} $n" >/dev/tty
        # Executa com LIKWID
        likwid-perfctr -O -C ${CPU} -g ${m} -o ${LIKWID_OUT} -m ./${PROG} ${n} >> ${TEMPOS}

        # Extrai métricas específicas para cada região (marcador LIKWID)
        matVet_value=$(grep -A 10 "Region matVet" ${LIKWID_OUT} | grep "${METRIC_PATTERN}" | awk '{print $NF}' | head -1)
        matVet_otim_value=$(grep -A 10 "Region matVet_otim" ${LIKWID_OUT} | grep "${METRIC_PATTERN}" | awk '{print $NF}' | head -1)
        matMat_value=$(grep -A 10 "Region matMat" ${LIKWID_OUT} | grep "${METRIC_PATTERN}" | awk '{print $NF}' | head -1)
        matmat_otim_value=$(grep -A 10 "Region matmat_otim" ${LIKWID_OUT} | grep "${METRIC_PATTERN}" | awk '{print $NF}' | head -1)
        
        # Trata valores vazios
        matVet_value=${matVet_value:-"0.0"}
        matVet_otim_value=${matVet_otim_value:-"0.0"}
        matMat_value=${matMat_value:-"0.0"}
        matmat_otim_value=${matmat_otim_value:-"0.0"}
        
        # Escreve linha no CSV
        echo "${n},${matVet_value},${matVet_otim_value},${matMat_value},${matmat_otim_value}" >> ${LIKWID_CSV}
    done
done

echo "powersave" > /sys/devices/system/cpu/cpufreq/policy${CPU}/scaling_governor 

