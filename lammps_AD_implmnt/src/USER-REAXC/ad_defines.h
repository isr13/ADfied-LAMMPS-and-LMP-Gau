#include <mpi.h>

#ifndef AD_DEFINES_H
#define AD_DEFINES_H
  #include "ad.hpp"  // use stce's OpenSource AD tool
  AD_ALLOW_EXPLICIT_TYPE_CAST_TO ( int )
  AD_ALLOW_EXPLICIT_TYPE_CAST_TO ( long )

  // first order scalar adjoint type
  using ADmode = ad::ga1s<double>;
  using ADtype = ADmode::type;
  using myScalar = ADtype;

  // fix for MPI_Allreduce with AD datatype
  inline int MY_MPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
    if (datatype == MPI_DOUBLE){
      return MPI_Allreduce(sendbuf, recvbuf, count*sizeof(myScalar)/sizeof(double), MPI_DOUBLE, op, comm);
    }else{
      return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    }
  }


#endif
