// ============================================================================
// to be added
//
// Author: Jiayi Huang and Ju Liu
// Date: Sept. 20 2022
// ============================================================================
#include "yaml.h"
#include "Vector_3.hpp"
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <cmath>

// normalize the vectors for comfortable visualization
int main(int argc, char *argv[])
{
  YAML::Node config = YAML::LoadFile("../config.yaml");
  std::string filename_vtu = config["solid_vtu"].as<std::string>();      // wall_vtu
  std::string filename_vtp = config["centerline_vtp"].as<std::string>(); // centerline_vtp

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

  // construct three gird containers to save three normals respectively
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
  vtkNew<vtkPointLocator> pt_locator;
  pt_locator->Initialize();
  pt_locator->SetDataSet(poly);
  pt_locator->BuildLocator();

  // Find the closest points
  double pt[3];
  double cl_pt_1[3];
  double cl_pt_2[3];

  // Initiate vectors for the computation of distance and normalization
  Vector_3 pt1, pt2, uu, vv, ww;
  int cl_pt_Id;

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
  double u[3], v[3];

  for (int ii = 0; ii < nFunc; ++ii)
  {

    grid->GetPoint(ii, pt);

    cl_pt_Id = pt_locator->FindClosestPoint(&pt[0]);
    polypoints->GetPoint(cl_pt_Id, cl_pt_1);

    // radial
    for (int jj = 0; jj < 3; ++jj)
    {
      u[jj] = pt[jj] - cl_pt_1[jj];
    }

    // longitudinal
    polypoints->GetPoint(cl_pt_Id - 1, cl_pt_1);
    polypoints->GetPoint(cl_pt_Id + 1, cl_pt_2);

    pt1.copy(cl_pt_1);
    pt2.copy(cl_pt_2);

    if (dist(pt1, pt2) >= 1)
    {
      polypoints->GetPoint(cl_pt_Id, cl_pt_1);
      polypoints->GetPoint(cl_pt_Id, cl_pt_2);
    }

    for (int jj = 0; jj < 3; ++jj)
    {
      v[jj] = cl_pt_1[jj] - cl_pt_2[jj];
    }

    uu.copy(u);
    vv.copy(v);

    // cross to get circumferential
    ww.copy(cross_product(uu, vv));

    uu.normalize();
    vv.normalize();
    ww.normalize();

    radial_normal->InsertTuple3(ii, uu.x(), uu.y(), uu.z());
    longitudinal_normal->InsertTuple3(ii, vv.x(), vv.y(), vv.z());
    circumferential_normal->InsertTuple3(ii, ww.x(), ww.y(), ww.z());
    if (ii % 10000 == 0)
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

  return 0;
}

// EOF
