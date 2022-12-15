// ============================================================================
// to be added
//
// Author: Jiayi Huang and Ju Liu
// Date: Sept. 20 2022
// ============================================================================
#include "yaml-cpp/yaml.h"
#include "Vector_3.hpp"
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkKdTreePointLocator.h>
#include <vtkOBBTree.h>
#include <cmath>

int main(int argc, char *argv[])
{
  YAML::Node config = YAML::LoadFile("../config.yaml");
  std::string filename_vtu = config["solid_vtu"].as<std::string>();                              // wall_vtu need to determine the fiber directions
  std::string filename_vtp = config["centerline_vtp"].as<std::string>();                         // centerline_vtp which has been modified(main cut branch)
  std::string lumen_wall_vtp = config["lumen_wall_vtp"].as<std::string>();                       // lumen_wall_vtp used for normal directions
  std::string tissue_wall_vtp = config["tissue_wall_vtp"].as<std::string>();                     // tissue_wall_vtp used for normal directions
  std::string local_radius_scale_factor = config["local_radius_scale_factor"].as<std::string>(); // to accelerate the search process
  std::string dot_product = config["dot_product"].as<std::string>();                             // adjust to find more correct centerline point

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
  int npoly = poly->GetNumberOfPoints();
  vtkPoints *polypoints = poly->GetPoints();
  vtkDataArray *radiusArray = poly->GetPointData()->GetArray("MaximumInscribedSphereRadius");

  vtkNew<vtkXMLPolyDataReader> reader_lumen;
  reader_lumen->SetFileName(lumen_wall_vtp.c_str());
  reader_lumen->Update();
  vtkPolyData *lumen = reader_lumen->GetOutput();
  vtkPoints *lumenpoints = lumen->GetPoints();
  vtkDataArray *normal_lumen = lumen->GetPointData()->GetArray("Normals");

  vtkNew<vtkXMLPolyDataReader> reader_tissue;
  reader_tissue->SetFileName(tissue_wall_vtp.c_str());
  reader_tissue->Update();
  vtkPolyData *tissue = reader_tissue->GetOutput();
  vtkPoints *tissuepoints = tissue->GetPoints();
  vtkDataArray *normal_tissue = tissue->GetPointData()->GetArray("Normals");

  // Identify the closest point on the centerline
  vtkNew<vtkPointLocator> pt_locator;
  pt_locator->Initialize();
  pt_locator->SetDataSet(poly);
  pt_locator->BuildLocator();

  vtkNew<vtkKdTreePointLocator> pt_radius_locator;
  pt_radius_locator->Initialize();
  pt_radius_locator->SetDataSet(poly);
  pt_radius_locator->BuildLocator();

  vtkNew<vtkPointLocator> pt_lumen_locator;
  pt_lumen_locator->Initialize();
  pt_lumen_locator->SetDataSet(lumen);
  pt_lumen_locator->BuildLocator();

  vtkNew<vtkPointLocator> pt_tissue_locator;
  pt_tissue_locator->Initialize();
  pt_tissue_locator->SetDataSet(tissue);
  pt_tissue_locator->BuildLocator();

  double pt[3];
  double ppt[3];
  double cl_pt_1[3];
  double cl_pt_2[3];

  // Initiate vectors for the computation of distance and normalization
  Vector_3 pt1, pt2, pt3, uu, vv, ww;
  int cl_pt_Id;

  // Three normal vectors definition
  // construct three gird containers to save three normals respectively
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

  // Create the locator
  vtkNew<vtkOBBTree> tree;
  tree->SetDataSet(reader_lumen->GetOutput());
  tree->BuildLocator();

  vtkNew<vtkOBBTree> tree_tissue;
  tree_tissue->SetDataSet(reader_tissue->GetOutput());
  tree_tissue->BuildLocator();

  double u[3], v[3], out1[3], out2[3];
  int cnt = 0, lumen_cl, tissue_cl;
  double tol = 0.1;

  double tmp_1, tmp_2;
  double ratio_lu, ratio_ti;
  vtkNew<vtkIdList> centerline_search;
  int number_radius = 0;
  int point_ind = 0;
  double input_radius = 0;

  for (int ii = 0; ii < nFunc; ii++)
  {
    cl_pt_Id = -1;
    grid->GetPoint(ii, pt);

    // radial
    lumen_cl = pt_lumen_locator->FindClosestPoint(pt);
    tissue_cl = pt_tissue_locator->FindClosestPoint(pt);

    lumenpoints->GetPoint(lumen_cl, out1);
    tissuepoints->GetPoint(tissue_cl, out2);

    for (int jj = 0; jj < 3; ++jj)
    {
      u[jj] = pt[jj] - out1[jj];
      v[jj] = out2[jj] - pt[jj];
    }

    pt1.copy(u);
    pt2.copy(v);

    for (int jj = 0; jj < 3; ++jj)
    {
      u[jj] = *(normal_lumen->GetTuple3(lumen_cl) + jj);
      v[jj] = *(normal_tissue->GetTuple3(tissue_cl) + jj);
    }

    uu.copy(u);
    vv.copy(v);
    uu.normalize();
    vv.normalize();

    tmp_1 = pt1.dot_product(uu) * pt1.dot_product(uu);
    tmp_2 = pt2.dot_product(vv) * pt2.dot_product(vv);
    ratio_lu = tmp_2 / (tmp_1 + tmp_2);
    ratio_ti = tmp_1 / (tmp_1 + tmp_2);
    uu.scale(ratio_lu);
    vv.scale(ratio_ti);
    if (uu.dot_product(vv) < 0)
      uu = vv - uu;
    else
      uu += vv;
    uu.normalize();

    radial_normal->InsertTuple3(ii, uu.x(), uu.y(), uu.z());

    // longitudinal
    poly->GetPoint(pt_locator->FindClosestPoint(pt), ppt);
    pt1.copy(pt);
    pt2.copy(ppt);
    input_radius = stod(local_radius_scale_factor) * dist(pt1, pt2); // determine the search region radius
    pt_radius_locator->FindPointsWithinRadius(input_radius, pt, centerline_search); // the Id of points within radius

    tmp_1 = 1000;
    tmp_2 = -1000;

    number_radius = centerline_search->GetNumberOfIds();

    for (int jj = 0; jj < number_radius; ++jj)
    {
      point_ind = centerline_search->GetId(jj);
      poly->GetPoint(point_ind, ppt);
      pt1.copy(pt);
      pt2.copy(ppt);
      pt3 = pt1 - pt2;
      pt3.normalize();

      if (uu.dot_product(pt3) > stod(dot_product))
      {
        if ((dist(pt1, pt2) / radiusArray->GetComponent(point_ind, 0)) < tmp_1)
        {
          cl_pt_Id = point_ind;
          tmp_1 = dist(pt1, pt2) / radiusArray->GetComponent(point_ind, 0);
        }
      }
    }

    if (cl_pt_Id > npoly - 1 || cl_pt_Id < 0)
    {
      printf("failed for no closetest point\n");
    }

    polypoints->GetPoint(cl_pt_Id + 1, cl_pt_1);
    polypoints->GetPoint(cl_pt_Id - 1, cl_pt_2);

    pt1.copy(cl_pt_1);
    pt2.copy(cl_pt_2);

    if (dist(pt1, pt2) >= tol)
    {
      polypoints->GetPoint(cl_pt_Id + 1, cl_pt_1);
      polypoints->GetPoint(cl_pt_Id, cl_pt_2);
    }

    pt1.copy(cl_pt_1);
    pt2.copy(cl_pt_2);

    if (dist(pt1, pt2) >= tol)
    {
      polypoints->GetPoint(cl_pt_Id, cl_pt_1);
      polypoints->GetPoint(cl_pt_Id - 1, cl_pt_2);
    }

    pt1.copy(cl_pt_1);
    pt2.copy(cl_pt_2);

    if (dist(pt1, pt2) >= tol)
    {
      printf("failed for too far\n");
      return 0;
    }

    // not correct longitudinal
    vv.copy(pt1 - pt2);
    
    // circumferential
    ww.copy(cross_product(uu, vv));
    ww.normalize();

    // correct longitudinal
    vv.copy(cross_product(uu, ww));
    vv.normalize();

    longitudinal_normal->InsertTuple3(ii, vv.x(), vv.y(), vv.z());
    circumferential_normal->InsertTuple3(ii, ww.x(), ww.y(), ww.z());
    if (ii % 1000 == 0)
      printf("normal ii = %d\n", ii);
  }
  printf("total undifine points = %d\n", cnt);
  grid->GetPointData()->AddArray(radial_normal);
  grid->GetPointData()->AddArray(longitudinal_normal);
  grid->GetPointData()->AddArray(circumferential_normal);

  vtkNew<vtkXMLUnstructuredGridWriter> writer;
  std::string str_before = "out";
  std::string str_after = filename_vtu.c_str();
  std::string str_combine = str_before + str_after;
  writer->SetFileName(str_combine.c_str());
  writer->SetInputData(grid);
  writer->Write();

  return 0;
}

// EOF