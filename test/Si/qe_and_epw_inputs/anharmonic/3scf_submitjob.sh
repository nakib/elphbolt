#!/bin/bash

INPDIR=$(pwd)/anharm_sc555_6nn
cd $INPDIR

for a in $(seq -f %03.0f 1 1 164) 
do
  echo $a
  SCFINPn="DISP.si_sc.in.${a}"
  n="DISP_3RD.${a}"
  echo $n
  echo $SCFINPn

  mkdir $n
  cd $n

  cp ../$SCFINPn ./scf.in

cat>Run.pbs<<EOF
#!/bin/bash
#PBS -l pmem=400mb,nodes=1:ppn=4,walltime=1:00:00:00
#PBS -N si.${a}

source /usr/public/Modules/default/init/bash
source /usr/public/Modules/default/init/bash_completion

cd \$PBS_O_WORKDIR
module load intel/2013.1.117 openmpi/1.8.4.intel

mpirun -np 4 pw.x < scf.in > scf.out

EOF

chmod 750 Run.pbs
qsub Run.pbs

cd ..

done
