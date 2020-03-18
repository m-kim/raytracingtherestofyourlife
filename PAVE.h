#ifndef PAVE_H
#define PAVE_H
#include <string>
#include <sstream>

#include <adios2.h>
#include <vtkm/cont/ArrayHandle.h>
#include <mpi.h>

class PAVE
{
public:
  PAVE(std::string _fn):
    fn(_fn)
  {
    adios = adios2::ADIOS(MPI_COMM_WORLD);
    bpIO = adios.DeclareIO("BPFile_N2N");
    writer = bpIO.Open(fn, adios2::Mode::Write);
  }
  ~PAVE()
  {
    writer.Close();
  }


  void save(std::string vn,
                 int nx, int ny,
                 const vtkm::cont::ArrayHandle<vtkm::Float32> &cols)
  {
    using ArrayType = vtkm::Float32;

    adios2::Variable<ArrayType> bpOut = bpIO.DefineVariable<ArrayType>(
          vn, {}, {}, {static_cast<std::size_t>(nx*ny)}, adios2::ConstantDims);



    auto *ptr = cols.GetStorage().GetArray();
    writer.Put<vtkm::Float32>(bpOut, ptr, adios2::Mode::Sync );
  }

  void save(std::string vn,
                 int nx, int ny,
                 const vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,4>> &cols)
  {
    using OutArrayType = vtkm::Float32;

    adios2::Variable<OutArrayType> bpOut = bpIO.DefineVariable<OutArrayType>(
          vn, {}, {}, {static_cast<std::size_t>(nx*ny*4)}, adios2::ConstantDims);


    std::vector<OutArrayType> arrayOut(cols.GetNumberOfValues()*4);
    for (int i=0; i<cols.GetNumberOfValues(); i++){
      auto col = cols.ReadPortal().Get(i);
      arrayOut[i*4] = col[0];
      arrayOut[i*4 + 1] = col[1];
      arrayOut[i*4 + 2] = col[2];
      arrayOut[i*4 + 3] = col[3];
    }
    writer.Put<OutArrayType>(bpOut, arrayOut.data());
  }

  void save(std::string vn,
                 int nx, int ny,
                 vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Int8,4>> &cols)
  {
    using OutArrayType = vtkm::Int8;

    adios2::Variable<OutArrayType> bpOut = bpIO.DefineVariable<OutArrayType>(
          vn, {}, {}, {static_cast<std::size_t>(nx*ny*4)}, adios2::ConstantDims);


    std::vector<OutArrayType> arrayOut(cols.GetNumberOfValues()*4);
    for (int i=0; i<cols.GetNumberOfValues(); i++){
      auto col = cols.ReadPortal().Get(i);
      arrayOut[i*4] = col[0];
      arrayOut[i*4 + 1] = col[1];
      arrayOut[i*4 + 2] = col[2];
      arrayOut[i*4 + 3] = col[3];
    }
    writer.Put<OutArrayType>(bpOut, arrayOut.data(), adios2::Mode::Sync );
  }

  void flush()
  {
    writer.PerformPuts();
    writer.Flush();
  }

  void close()
  {
    writer.Close();
  }
  adios2::ADIOS adios;
  adios2::IO bpIO;
  adios2::Engine writer;

  std::string fn;

};

#endif
