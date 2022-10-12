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
#include <vtkKdTreePointLocator.h>
#include <cmath>

//normalize the vectors for comfortable visualization
int normalization(double *normal)
{
  double total = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
  for (int i = 0; i < 3; ++i)
    normal[i] /= total;
  return 0;
}

int main(int argc, char *argv[])
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  // Parse command line arguments: read in the wall_vtu and centerline_vtp to output local_ref_vtu
  if (argc != 3)
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
  vtkNew<vtkXMLUnstructuredGridReader> reader_vtu1;
  reader_vtu1->SetFileName(filename_vtu.c_str());
  reader_vtu1->Update();
  vtkNew<vtkXMLUnstructuredGridReader> reader_vtu2;
  reader_vtu2->SetFileName(filename_vtu.c_str());
  reader_vtu2->Update();

  //construct three gird containers to save three normals respectively
  vtkUnstructuredGrid *grid = reader_vtu->GetOutput();
  vtkUnstructuredGrid *grid1 = reader_vtu1->GetOutput();
  vtkUnstructuredGrid *grid2 = reader_vtu2->GetOutput();

  int nFunc = grid->GetNumberOfPoints();

  vtkNew<vtkXMLPolyDataReader> reader_vtp;
  reader_vtp->SetFileName(filename_vtp.c_str());
  reader_vtp->Update();
  vtkPolyData *poly = reader_vtp->GetOutput();
  vtkPoints *polypoints = poly->GetPoints();


  // Identify the closest point on the centerline
  vtkNew<vtkCellLocator> locator;
  locator->Initialize();
  locator->SetDataSet(poly);
  locator->BuildLocator();

  // Find the closest points
  vtkIdType cellId; // the cell id of the cell containing the closest point will be returned here
  double cl_pt[3];
  double nb_pt[3];
  vtkNew<vtkGenericCell> cell;
  int subId;
  double dist;

  // Create the tree
  vtkSmartPointer<vtkKdTreePointLocator> pointTree =
      vtkSmartPointer<vtkKdTreePointLocator>::New();
  pointTree->SetDataSet(poly);
  pointTree->BuildLocator();

  // Find the k closest points to (0.5, 0.5, 0.5)
  unsigned int k = 2;
  vtkSmartPointer<vtkIdList> result =
      vtkSmartPointer<vtkIdList>::New();

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
    locator->FindClosestPoint(&pt[0], &cl_pt[0], cell, cellId, subId, dist);//the clostest point
    pointTree->FindClosestNPoints(k, &pt[0], result);          //two clostest point 

    //radial
    for (int jj = 0; jj < 3; ++jj)
    {
      u[jj] = pt[jj] - cl_pt[jj];
      printf("pt[jj] - cl_pt[jj] = %lf - %lf\n", pt[jj], cl_pt[jj]);
    }
    
    //longitudinal
    polypoints->GetPoint(result->GetId(1), nb_pt);
    polypoints->GetPoint(result->GetId(0), cl_pt);

    for (int jj = 0; jj < 3; ++jj)
    {
      v[jj] = nb_pt[jj] - cl_pt[jj];
      printf("nb_pt[jj] - cl_pt[jj] = %lf - %lf\n", nb_pt[jj], cl_pt[jj]);
    }

    //cross to get circumferential
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];

    normalization(u);
    normalization(v);
    normalization(w);
    
    radial_normal->InsertTuple(ii, u);
    longitudinal_normal->InsertTuple(ii, v);
    circumferential_normal->InsertTuple(ii, w);
    printf("normal ii = %d\n", ii);
  }

  grid->GetPointData()->SetVectors(radial_normal);
  grid1->GetPointData()->SetVectors(longitudinal_normal);
  grid2->GetPointData()->SetVectors(circumferential_normal);

  vtkNew<vtkXMLUnstructuredGridWriter> writer;
  vtkNew<vtkXMLUnstructuredGridWriter> writer1;
  vtkNew<vtkXMLUnstructuredGridWriter> writer2;

  writer->SetFileName("radial_normal.vtu");
  writer->SetInputData(grid);
  writer->Write();
  writer1->SetFileName("longitudinal_normal.vtu");
  writer1->SetInputData(grid1);
  writer1->Write();
  writer2->SetFileName("circumferential_normal.vtu");
  writer2->SetInputData(grid2);
  writer2->Write();

  PetscFinalize();
  return 0;
}

// EOF
