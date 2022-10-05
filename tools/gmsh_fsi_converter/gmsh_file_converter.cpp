// ============================================================================
// to be added
//
// Author: Jiayi Huang and Ju Liu
// Date: Sept. 20 2022
// ============================================================================
#include "Tet_Tools.hpp"
#include "ElemBC_3D_tet_wall.hpp"
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <fstream>

int main(int argc, char *argv[])
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  // Parse command line arguments: read in the wall_vtu and centerline_vtp to output local_ref_vtu
  if (argc != 4)
  {
    std::cerr << "Usage: " << argv[0] << " Filename to read wall_vtu"
              << " Filename to read centerline_vtp"
              << " Filename to write local_ref_vtu"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string filename_vtu = argv[1]; // wall_vtu
  std::string filename_vtp = argv[2]; // centerline_vtp

  // Read all the data from the file
  vtkNew<vtkXMLUnstructuredGridReader> reader_vtu;
  reader_vtu->SetFileName(filename_vtu.c_str());
  reader_vtu->Update();
  vtkUnstructuredGrid *grid = reader_vtu->GetOutput();
  int nFunc = grid->GetNumberOfPoints();

  vtkNew<vtkXMLPolyDataReader> reader_vtp;
  reader_vtp->SetFileName(filename_vtp.c_str());
  reader_vtp->Update();
  vtkPolyData *poly = reader_vtp->GetOutput();
  // int num_centerline = poly->GetNumberOfPoints();

  // Identify the closest point on the centerline
  vtkCellLocator *locator = vtkCellLocator::New();
  locator->Initialize();
  locator->SetDataSet(poly);
  locator->BuildLocator();

  // Find the closest points
  vtkIdType cellId; // the cell id of the cell containing the closest point will be returned here
  double cl_pt[3];
  double nb_pt[3];
  vtkGenericCell *cell = vtkGenericCell::New();
  int subId;
  double dist;

  // Three normal vectors definition
  vtkNew<vtkDoubleArray> radial_normal, longitudinal_normal, circumferential_normal;
  radial_normal->SetName("radial_normal");
  radial_normal->SetNumberOfComponents(3);
  radial_normal->SetNumberOfTuples(nFunc);
  longitudinal_normal->SetName("longitudinal_normal");
  longitudinal_normal->SetNumberOfComponents(3);
  longitudinal_normal->SetNumberOfTuples(nFunc);
  circumferential_normal->SetName("circumferential_normal");
  circumferential_normal->SetNumberOfComponents(3);
  circumferential_normal->SetNumberOfTuples(nFunc);
  double u[3], v[3], w[3];

  for (int ii = 0; ii < nFunc; ++ii)
  {
    double pt[3];
    grid->GetPoint(ii, pt);
    locator->FindClosestPoint(&pt[0], &cl_pt[0], cell, cellId, subId, dist);
    locator->FindClosestPoint(&cl_pt[0], &nb_pt[0], cell, cellId, subId, dist);
    for (int jj = 0; jj < 3; ++jj)
    {
      u[jj] = pt[jj] - cl_pt[jj];
      v[jj] = cl_pt[jj] - nb_pt[jj];
    }
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
    radial_normal->InsertTuple(ii, u);
    printf("radial_normal ii = %d\n", ii);
    longitudinal_normal->InsertTuple(ii, v);
    printf("longitudinal_normal ii = %d\n", ii);
    circumferential_normal->InsertTuple(ii, w);
    printf("circumferential_normal ii = %d\n", ii);
  }

  grid->GetPointData()->SetVectors(radial_normal);
  // grid->GetPointData()->SetVectors(longitudinal_normal);
  // grid->GetPointData()->SetVectors(circumferential_normal);
  std::string write_vtu = argv[3];
  vtkNew<vtkXMLUnstructuredGridWriter> writer;
  writer->SetFileName(write_vtu.c_str());
  writer->SetInputData(grid);
  writer->Write();

  // clean memory
  locator->Delete();
  cell->Delete();

  PetscFinalize();
  return 0;
}

// EOF
