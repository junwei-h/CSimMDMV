#!/bin/csh 

#nonconditional simulation, grid interval 2 m
set num=$1
while ( $num >= $3 )
  mpirun -np $num ../bin/Sim3D3V_BENCH  ./bench/grids200/Sim3D3V_ds200.inp >./bench/grids200/Sim3D3V_ds200_np$num.out
  echo "Sim3D3V_ds200_np$num is done"
  @ num = $num - $2
end

#nonconditional simulation, grid interval 1 m
set num=$1
while ( $num >= $3 )
  mpirun -np $num ../bin/Sim3D3V_BENCH  ./bench/grids100/Sim3D3V_ds100.inp >./bench/grids100/Sim3D3V_ds100_np$num.out
  echo "Sim3D3V_ds100_np$num is done"
  @ num = $num - $2
end

#conditional simulation, grid interval 2 m
set num=$1
while ( $num >= $3 )
  mpirun -np $num ../bin/Sim3D3V_BENCH  ./bench/grids200/CSim3D3V_ds200.inp >./bench/grids200/CSim3D3V_ds200_np$num.out
  echo "CSim3D3V_ds200_np$num is done"
  @ num = $num - $2
end

#conditional simulation, grid interval 1 m
set num=$1
while ( $num >= $3 )
  mpirun -np $num ../bin/Sim3D3V_BENCH  ./bench/grids100/CSim3D3V_ds100.inp >./bench/grids100/CSim3D3V_ds100_np$num.out
  echo "CSim3D3V_ds100_np$num is done"
  @ num = $num - $2
end

#nonconditional simulation, grid interval 0.4 m
set num=$1
while ( $num >= $3 )
  mpirun -np $num ../bin/Sim3D3V_BENCH  ./bench/grids40/Sim3D3V_ds40.inp >./bench/grids40/Sim3D3V_ds40_np$num.out
  echo "Sim3D3V_ds40_np$num is done"
  @ num = $num - $2
end

#conditional simulation, grid interval 0.4 m
set num=$1
while ( $num >= $3 )
  mpirun -np $num ../bin/Sim3D3V_BENCH  ./bench/grids40/CSim3D3V_ds40.inp >./bench/grids40/CSim3D3V_ds40_np$num.out
  echo "CSim3D3V_ds40_np$num is done"
  @ num = $num - $2
end