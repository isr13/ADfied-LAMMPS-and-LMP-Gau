#include <mpi.h>

#ifndef AD_DEFINES_H
#define AD_DEFINES_H
  #include "ad.hpp"
  AD_ENABLE_TYPE_CONSTRUCTION_FROM ( int )
  AD_ENABLE_TYPE_CONSTRUCTION_FROM ( long )
  typedef ad::ga1s <double >::type myScalar ;
  //typedef double myScalar;
  //typedef ad::gt1s<double>::type myScalar;
  
  inline int MY_MPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
    if (datatype == MPI_DOUBLE){
      return MPI_Allreduce(sendbuf, recvbuf, count*sizeof(myScalar)/sizeof(double), MPI_DOUBLE, op, comm);
    }else{
      return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    }
  }


#endif
