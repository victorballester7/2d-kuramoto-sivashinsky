#include <getopt.h>

#include "mpitranspose.h"
#include "mpiutils.h"
#include "statistics.h"
#include "align.h"
#include "timing.h"

using namespace std;
using namespace utils;

//namespace utils {
//unsigned int defaultmpithreads=1;
//}

const size_t showlimit=1024;
size_t N=1000000;

void init(Complex *data, size_t X, size_t y, size_t Z,
          int x0, int y0) {
  for(size_t i=0; i < X; ++i) {
    for(size_t j=0; j < y; ++j) {
      for(size_t k=0; k < Z; ++k) {
        data[(y*i+j)*Z+k].re=x0+i;
        data[(y*i+j)*Z+k].im=y0+j;
      }
    }
  }
}

inline void usage()
{
  cerr << "Options: " << endl;
  cerr << "-h\t\t help" << endl;
  cerr << "-T<int>\t\t number of threads" << endl;
  cerr << "-t\t\t test" << endl;
  cerr << "-N<int>\t\t number of timing tests" << endl;
  cerr << "-m<int>\t\t size" << endl;
  cerr << "-x<int>\t\t x size" << endl;
  cerr << "-y<int>\t\t y size" << endl;
  cerr << "-z<int>\t\t z size" << endl;
  cerr << "-S<int>\t\t stats choice" << endl;
  cerr << "-p<int>\t\t which part of the transpose to time" << endl;
  usageTranspose();
  cerr << "-L\t\t locally transpose output" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  bool test=false;
  bool quiet=false;

  size_t X=8, Y=8, Z=1;
  int a=0; // Test for best block divisor
  int alltoall=-1; // Test for best alltoall routine


  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);

  int stats=0;

  int size,rank;

  MPI_Comm communicator=MPI_COMM_WORLD;
  MPI_Comm_size(communicator,&size);
  MPI_Comm_rank(communicator,&rank);

  bool main=rank == 0;

  if(!main)
    opterr=0;

#ifdef __GNUC__
  optind=0;
#endif
  for (;;) {
    int c=getopt(argc,argv,"h:A:a:m:N:s:T:S:x:y:z:qt");
    if (c == -1) break;

    switch (c) {
      case 0:
        break;
      case 's':
        alltoall=atoi(optarg);
        break;
      case 'a':
        a=atoi(optarg);
        break;
      case 'm':
        X=Y=atoi(optarg);
        break;
      case 'x':
        X=atoi(optarg);
        break;
      case 'y':
        Y=atoi(optarg);
        break;
      case 'z':
        Z=atoi(optarg);
        break;
      case 'T':
        defaultmpithreads=atoi(optarg);
        break;
      case 'S':
        stats=atoi(optarg);
        break;
      case 't':
        test=true;
        break;
      case 'q':
        quiet=true;
        break;
      case 'N':
        N=atoi(optarg);
        break;
      case 'h':
      default:
        if(rank == 0)
          usage();
    }
  }

  if(provided < MPI_THREAD_FUNNELED)
    defaultmpithreads=1;

  if(rank == 0) {
    cout << "size=" << size << endl;
    cout << "threads=" << defaultmpithreads << endl;
  }

  if(test) N=1;

  int retval=0;

  split d(X,Y,communicator);

  size_t x=d.x;
  size_t y=d.y;
  size_t y0=d.y0;

  if(main) {
    cout << "size=" << size << endl;
    cout << "x=" << x << endl;
    cout << "y=" << y << endl;
    cout << "X=" << X << endl;
    cout << "Y=" << Y << endl;
    cout << "Z=" << Z << endl;
    cout << "N=" << N << endl;
  }

  Complex* data=ComplexAlign(max(X*y,x*Y)*Z);

  init(data,X,y,Z,0,d.y0);

  //    show(data,X,y*Z,communicator);

  mpitranspose<Complex> T(X,Y,x,y,Z,data,NULL,communicator,
                          mpiOptions(a,alltoall,defaultmpithreads,!quiet));
  init(data,X,y,Z,0,y0);
  T.localize1(data);

  //    show(data,x,Y*Z,communicator);

  init(data,X,y,Z,0,y0);

  statistics Sininit,Sinwait0,Sinwait1,Sin;
  statistics Soutinit,Soutwait0,Soutwait1,Sout;

  bool showoutput=!quiet && (test || (!X*Y < showlimit && N == 1));

  if(showoutput) {
    if(main)
      cout << "\nInput:" << endl;
    show(data,X,y*Z,communicator);
  }


  if(test) {
    if(main)
      cout << "\nDiagnostics and unit test.\n" << endl;

    init(data,X,y,Z,0,y0);
    if(showoutput) {
      if(main)
        cout << "Input:" << endl;
      show(data,X,y*Z,communicator);
    }

    Complex *wholedata=NULL, *wholeoutput=NULL;
    if(main) {
      wholedata=new Complex[X*Y*Z];
      wholeoutput=new Complex[X*Y*Z];
    }

    gathery(data,wholedata,d,Z,communicator);

    if(showoutput && main) {
      cout << "\nGathered input data:" << endl;
      show(wholedata,X,Y,0,0,X,Y);
    }

    T.localize1(data); // N x m -> n x M
    T.localize0(data); // n x M -> N x m
    T.localize1(data); // N x m -> n x M

    if(showoutput) {
      if(main)
        cout << "\nOutput:" << endl;
      show(data,X,y*Z,communicator);
    }

    gatherx(data,wholeoutput,d,Z,communicator);

    if(main) {
      if(showoutput) {
        cout << "\nGathered output data:" << endl;
        show(wholeoutput,X,Y,0,0,X,Y);
      }

      bool success=true;
      const size_t stop=X*Y*Z;
      for(size_t pos=0; pos < stop; ++pos) {
        if(wholedata[pos] != wholeoutput[pos])
          success=false;
      }

      if(success == true) {
        cout << "\nTest succeeded." << endl;
      } else {
        cout << "\nERROR: TEST FAILED!" << endl;
        ++retval;
      }

    }
  } else {
    if(main)
      cout << "\nSpeed test.\n" << endl;

    bool detailed=main && stats == -1;
    vector<double> Tp;
    init(data,X,y,Z,0,y0);
    T.localize0(data); // Initialize communication buffers

    for(size_t k=0; k < N; ++k) {
      init(data,X,y,Z,0,y0);

      double begin=0.0, Tinit0=0.0, Tinit=0.0, Twait0=0.0, Twait1=0.0;
      utils::stopWatch *W=NULL;
      if(main) W=new stopWatch;
      T.inphase0();
      if(main) Tinit0=W->seconds();
      T.insync0();
      if(main) Twait0=W->seconds();
      T.inphase1();
      if(main) Tinit=W->seconds();
      T.insync1();
      if(main) Twait1=W->seconds();
      T.inpost();

      double tin=0.0;
      if(main) {
        tin=W->seconds();
        Sin.add(tin);
        Sininit.add(Tinit0-begin);
        Sinwait0.add(Twait0-Tinit0);
        Sinwait1.add(Twait1-Tinit);
      }

      if(showoutput) {
        if(main) cout << "Transpose:" << endl;
        show(data,x,Y*Z,communicator);
        if(main) cout << endl;
      }

      if(main) begin=W->seconds();
      T.outphase0();
      if(main) Tinit0=W->seconds();
      T.outsync0();
      if(main) Twait0=W->seconds();
      T.outphase1();
      if(main) Tinit=W->seconds();
      T.outsync1();
      if(main) Twait1=W->seconds();

      if(main) {
        double tout=W->seconds()-begin;
        if(detailed)
          Tp.push_back(0.5*(tout+tin));

        Sout.add(tout);
        Soutinit.add(Tinit0-begin);
        Soutwait0.add(Twait0-Tinit0);
        Soutwait1.add(Twait1-Tinit);
      }
    }

    if(detailed) {
      timings("transpose",X,Tp.data(),Tp.size(),stats);
      Tp.clear();
    }

    if(showoutput) {
      if(main) cout << "\nOriginal:" << endl;
      show(data,X,y*Z,communicator);
    }

    if(main) {
      cout << endl;
      Sininit.output("Tininit",X);
      Sinwait0.output("Tinwait0",X);
      Sinwait1.output("Tinwait1",X);
      Sin.output("Tin",X);
      cout << endl;
      Soutinit.output("Toutinit",X);
      Soutwait0.output("Toutwait0",X);
      Soutwait1.output("Toutwait1",X);
      Sout.output("Tout",X);
    }
  }

  deleteAlign(data);

  MPI_Finalize();
  return retval;
}
