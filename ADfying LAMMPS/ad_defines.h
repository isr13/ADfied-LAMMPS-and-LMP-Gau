#include <mpi.h>

#ifndef AD_DEFINES_H
#define AD_DEFINES_H
  #define DCO_DISABLE_AVX2_WARNING
  #define DCO_AUTO_SUPPORT
  #include "/home/em181266/lammps_dco/dco_cpp/dcl6i34ngl/include/dco.hpp"  // change for dco path
  DCO_ENABLE_EXPLICIT_TYPE_CAST_TO ( int )
  DCO_ENABLE_EXPLICIT_TYPE_CAST_TO ( long )
  typedef dco :: ga1s <double >:: type myScalar ;
  //typedef double myScalar;
  //typedef dco::gt1s<double>::type myScalar;
  
  inline int MY_MPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
    if (datatype == MPI_DOUBLE){
      return MPI_Allreduce(sendbuf, recvbuf, count*sizeof(myScalar)/sizeof(double), MPI_DOUBLE, op, comm);
    }else{
      return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    }
  }


#endif
