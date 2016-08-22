#!/bin/bash

set -e
set -u

name="$1"
overlapped=${name}-membrane-overlap
system=${name}-membrane
mdpPath=default_conf/

function convertToGMX {
    echo 1 3 | gmx pdb2gmx -f ${name}.pdb -o ${name}.gro -i ${name}-posre.itp -p ${name}.top
}

function alignProteinAndMembrane {
    lambada -f1 ${name}-gromacs-dimer.gro -f2 ../membrane/popc-1152-eq.gro -o ${overlapped}.gro
    rm rotated.*.gro ang*dat *oriented*
}

function removeOverlap {
    gmx editconf -f ${overlapped}.gro -o ${overlapped}.pdb
    pymol -d "
        load ${overlapped}.pdb;
        create protein, pol;
        create water, sol & (not (byres ((sol) nto. 1.5 of pol)))
        create popc, resn pop & (not (byres ((resn pop) nto. 0.75 of pol)))
        save ${system}.pdb, protein | water | popc;
        quit;
    "
}

function updateTop {
    npopc=`grep C52 ${system}.pdb | wc -l`
    nsol=`grep HW1 ${system}.pdb | wc -l`

    top="${system}.top"
    sed '/\[\s*molecules\s*\]/,$d' ${name}.top > ${top}
    echo '#include "./oplsaa.ff/popc_opls.itp"' >> ${top}
    echo "[ molecules ]" >> ${top}
    echo "Protein_chain_A 2" >> ${top}
    echo "SOL ${nsol}" >> ${top}
    echo "POPC ${npopc}" >> ${top}
}

function runMd {
    base=$1
    cpi=''
    if [ -f ${base}.cpt ] ; then
        cpi="-cpi ${base}.cpt"
    fi

    gmx mdrun -v -deffnm ${base} -cpt 5 ${cpi} -c ${base}.pdb
}

function runVacuumMin {
    gmx editconf -f ${system}.pdb -o ${system}.gro
    gmx grompp -f ${mdpPath}/min-vacuum.mdp -po min-vacuum.mdp -c ${system}.gro -p ${system}.top -o em-vacuum.tpr
    runMd em-vacuum
}

function solvate {
    gmx solvate -cp em-vacuum.pdb -o ${system}-solvated.pdb -p ${system}.top 
}

function runMin {
    gmx editconf -f ${system}-solvated.pdb -o ${system}-solvated.gro
    gmx grompp -f ${mdpPath}/min-vacuum.mdp -po min-vacuum.mdp -c ${system}-solvated.gro -p ${system}.top -o em.tpr
    runMd em
}

function runNVT {
    gmx grompp -f ${mdpPath}/nvt.mdp -po nvt.mdp  -c em.pdb -p ${system}.top -o nvt.tpr
    runMd nvt
}

function runNPT {
    gmx grompp -f ${mdpPath}/npt.mdp -po npt.mdp  -c nvt.pdb -p ${system}.top -o npt.tpr
    runMd npt
}

function runProd {
    gmx grompp -f ../default_conf/prod.mdp -po prod.mdp  -c npt.pdb -p ${system}.top -o prod.tpr
    runMd prod
}

alignProteinAndMembrane
removeOverlap
updateTop
runVacuumMin
solvate
runMin
runNVT
runNPT
runProd
