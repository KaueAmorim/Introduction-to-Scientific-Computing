#!/bin/bash

EXEC="./resolveEDO"
CSV="medidas.csv"
METRICA="FLOPS_DP"
NUCLEO=3

FREQUENCIA="/sys/devices/system/cpu/cpufreq/policy${NUCLEO}/scaling_governor"
if [ -w "${FREQUENCIA}" ]; then
    echo "performance" > ${FREQUENCIA}
fi

if [ -f ${EXEC} ]; then
    make purge
fi

make likwid

likwid-perfctr -C ${NUCLEO} -g ${METRICA} -O -o ${CSV} -m ${EXEC}
echo ""
grep FP_ARITH_INST_RETIRED_SCALAR_DOUBLE ${CSV} | cut -d ',' -f 1,3

if [ -f ${CSV} ]; then
    rm ${CSV}
fi

if [ -w "${FREQUENCIA}" ]; then
    echo "powersave" > ${FREQUENCIA}
fi