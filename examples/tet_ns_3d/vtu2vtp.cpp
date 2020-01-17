#include "vtkGeometryFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLUnstructuredGridReader.h"

#include <iostream>

int main(int argc, char* argv[])
{
  if (argc != 3)
  {
    std::cerr << "Usage:\nvtu2vtp input.vtu output.vtp\n";
    return 1;
  }

  // build unstructured grid reader
  auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

  // pass the filename along to the reader
  reader->SetFileName(argv[1]);

  // Force reading
  reader->Update();

  // read in the grid
  auto ugrid = reader->GetOutput();

  // process grid to poly
  auto geometryFilter = vtkGeometryFilter::New();
  geometryFilter->SetInputData(ugrid);
  geometryFilter->Update();
  auto polydata = geometryFilter->GetOutput();

  // set up writer
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(argv[2]);
  writer->SetInputData(polydata);
  writer->Write();

  return 0;
}

