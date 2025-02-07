// ============================================================================
// This is a file converter driver that reads in the vessel_wall vtu file format,
// the centerline vtp, the lumen_wall vtp, and the tissue_wall vtp files,
// and then define the local coordinate system of points on the vessel wall grid 
// according to the centerline.
// 
// This drive will read in the files listed in yaml.
// The yaml will have 4 lines.
//
// The following is very IMPORTANT!! 
// Make sure the first line is the vtu volume file.
// The second line represents the centerline vtp file, 
// created by merging multiple branch centerline files, 
// with each centerline having undergone trimming.
// The third line refers to the lumen_wall vtp file, 
// which contains the outward-facing surface normal vectors for each point.
// The fourth line refers to the tissue_wall vtp file, 
// which contains the outward-facing surface normal vectors for each point.
//
// Author: Jiayi Huang and Ju Liu
// Date: Sept. 20 2022
// ============================================================================

#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "VTK_Tools.hpp"
#include <vtkDoubleArray.h>
#include <vtkParametricFunctionSource.h>
#include <vtkParametricSpline.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <yaml-cpp/yaml.h>

int main(int argc, char *argv[])
{
  // ===== Yaml options ===== 
  bool is_loadYaml = true;
  std::string yaml_file("../config.yml");

  // ===== Yaml Arguments =====
  SYS_T::GetOptionBool(  "-is_loadYaml",     is_loadYaml);
  SYS_T::GetOptionString("-yaml_file",       yaml_file);

  if (is_loadYaml) SYS_T::InsertFileYAML( yaml_file,  false );

  // ===== Yaml vtu,vtp Arguments =====
  YAML::Node config = YAML::LoadFile(yaml_file.c_str());
  std::string solid_vtu_s = config["solid_vtu"].as<std::string>();               // solid_wall_vtu
  std::string centerline_vtp_s = config["centerline_vtp"].as<std::string>();     // centerline_vtp
  std::string lumen_wall_vtp_s = config["lumen_wall_vtp"].as<std::string>();     // lumen_wall_vtp
  std::string tissue_wall_vtp_s = config["tissue_wall_vtp"].as<std::string>();   // tissue_wall_vtp

  // ===== Read Command Line Arguments =====
  SYS_T::commPrint("===> Reading command line arguments... \n");

  SYS_T::GetOptionString(   "-solid_vtu_s",         solid_vtu_s);
  SYS_T::GetOptionString(   "-centerline_vtp_s",    centerline_vtp_s);
  SYS_T::GetOptionString(   "-lumen_wall_vtp_s",    lumen_wall_vtp_s);
  SYS_T::GetOptionString(   "-tissue_wall_vtp_s",   tissue_wall_vtp_s);
  
  // // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint(   "-solid_vtu_s",         solid_vtu_s);
  SYS_T::cmdPrint(   "-centerline_vtp_s",    centerline_vtp_s);
  SYS_T::cmdPrint(   "-lumen_wall_vtp_s",    lumen_wall_vtp_s);
  SYS_T::cmdPrint(   "-tissue_wall_vtp_s",   tissue_wall_vtp_s);

  // ===== Read solid mesh ===== 
  int nFunc_s, nElem_s;
  std::vector<int> vecIEN_s;
  std::vector<double> ctrlPts_s;
  VTK_T::read_vtu_grid(solid_vtu_s, nFunc_s, nElem_s, ctrlPts_s, vecIEN_s);
  vtkNew<vtkXMLUnstructuredGridReader> solid_vtu;
  solid_vtu->SetFileName(solid_vtu_s.c_str());
  solid_vtu->Update();
  vtkUnstructuredGrid *solid_grid = solid_vtu->GetOutput();

  // ===== Read the lumen wall vtp ===== 
  // This lumen vtp is assumed to consist of all normal pointing out.
  vtkNew<vtkXMLPolyDataReader> lumen_wall_vtp;
  lumen_wall_vtp->SetFileName(lumen_wall_vtp_s.c_str());
  lumen_wall_vtp->Update();
  vtkPolyData *lumen_wall_poly = lumen_wall_vtp->GetOutput();
  vtkPoints *lumen_wall_points = lumen_wall_poly->GetPoints();
  vtkDataArray *lumen_wall_normal = lumen_wall_poly->GetPointData()->GetArray("Normals");

  // ===== Read the tissue wall vtp ===== 
  // This tissue vtp is assumed to consist of all normal pointing out.
  vtkNew<vtkXMLPolyDataReader> tissue_wall_vtp;
  tissue_wall_vtp->SetFileName(tissue_wall_vtp_s.c_str());
  tissue_wall_vtp->Update();
  vtkPolyData *tissue_wall_poly = tissue_wall_vtp->GetOutput();
  vtkPoints *tissue_wall_points = tissue_wall_poly->GetPoints();
  vtkDataArray *tissue_wall_normal = tissue_wall_poly->GetPointData()->GetArray("Normals");
  
  // ===== Three normal vectors definition =====  
  // Construct three arrays to save three local direction vecters. 
  // Letters r, l, and c denotes the radial, longitudinal,
  // and circumferential directions, respectively.
  vtkNew<vtkDoubleArray> loc_r_vec_array, loc_l_vec_array, loc_c_vec_array;
  loc_r_vec_array->SetName("loc_r_vec");
  loc_r_vec_array->SetNumberOfComponents(3);
  loc_r_vec_array->SetNumberOfTuples(nFunc_s);
  loc_l_vec_array->SetName("loc_l_vec");
  loc_l_vec_array->SetNumberOfComponents(3);
  loc_l_vec_array->SetNumberOfTuples(nFunc_s);
  loc_c_vec_array->SetName("loc_c_vec");
  loc_c_vec_array->SetNumberOfComponents(3);
  loc_c_vec_array->SetNumberOfTuples(nFunc_s);

  // ===== Set the locator of lumen surface vtp for finding the closest point =====
  vtkNew<vtkPointLocator> lumen_locator;
  lumen_locator->Initialize();
  lumen_locator->SetDataSet(lumen_wall_poly);
  lumen_locator->BuildLocator();

  // ===== Set the locator of tissue surface vtp for finding the closest point =====
  vtkNew<vtkPointLocator> tissue_locator;
  tissue_locator->Initialize();
  tissue_locator->SetDataSet(tissue_wall_poly);
  tissue_locator->BuildLocator();


  // ===== Read the centerline vtp ===== 
  // This centerline vtp is assumed to consist of all branch vessel centerline vtp.
  vtkNew<vtkXMLPolyDataReader> centerline_vtp;
  centerline_vtp->SetFileName(centerline_vtp_s.c_str());
  centerline_vtp->Update();
  vtkPolyData *centerline_poly = centerline_vtp->GetOutput();
  vtkPointData* centerline_pointData = centerline_poly->GetPointData();

  // ===== Create two containers to store the branches vtp ===== 
  // PolyDatas contains the the points coordinates and Locators are used to find the ClosestPoint.  
  std::vector<vtkSmartPointer<vtkPolyData>> PolyDatas;
  std::vector<vtkSmartPointer<vtkPointLocator>> Locators;

  // ===== Create one container to store relevant distance =====   
  std::vector<double> dist_relevant;
  for (int ii=0; ii<nFunc_s; ++ii) dist_relevant.push_back(10);

  // ===== Traverse all data attributes and find all attributes starting with "radius_" =====
  for (int ii = 0; ii < centerline_pointData->GetNumberOfArrays(); ++ii) {
      std::string arrayName = centerline_pointData->GetArrayName(ii);
      if (arrayName.find("radius_") == 0) {
          vtkDataArray* dataArray = centerline_pointData->GetArray(arrayName.c_str());

          double range[2];
          centerline_pointData->GetRange(arrayName.c_str(),range,0);

          vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();

          for (vtkIdType jj = 0; jj < dataArray->GetNumberOfTuples(); ++jj)
          {
            double pt[3];
            centerline_poly->GetPoint(jj + range[0], pt);
            newPoints->InsertNextPoint(pt);
          }

          vtkSmartPointer<vtkPolyData> newPolyData = vtkSmartPointer<vtkPolyData>::New();
          newPolyData->SetPoints(newPoints);
          newPolyData->GetPointData()->AddArray(dataArray);

          vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
          locator->SetDataSet(newPolyData);
          locator->BuildLocator();

          PolyDatas.push_back(newPolyData);
          Locators.push_back(locator);
      }
  }

  for (int ii = 0; ii < nFunc_s; ++ii)
  {
    // ===== Read the ii_th point of solid mesh ===== 
    double solid_pt[3]; 
    solid_grid->GetPoint(ii, solid_pt);

    // ===== Define and get the closest point from solid point to two surfaces =====
    double lumen_closest_solid_pt[3], tissue_closest_solid_pt[3]; 
    int id_lumen_closest_solid_pt=lumen_locator->FindClosestPoint(solid_pt);
    int id_tissue_closest_solid_pt=tissue_locator->FindClosestPoint(solid_pt);    
    lumen_wall_points->GetPoint(id_lumen_closest_solid_pt, lumen_closest_solid_pt);
    tissue_wall_points->GetPoint(id_tissue_closest_solid_pt, tissue_closest_solid_pt);

    // ===== Use Vector_3 store the closest point =====
    Vector_3 lumen_closest_solid_vec(lumen_closest_solid_pt[0],lumen_closest_solid_pt[1],lumen_closest_solid_pt[2]);
    Vector_3 tissue_closest_solid_vec(tissue_closest_solid_pt[0],tissue_closest_solid_pt[1],tissue_closest_solid_pt[2]);
    Vector_3 solid_vec(solid_pt[0],solid_pt[1],solid_pt[2]);

    // ===== Compute the distance from solid point to the closest point =====
    double dist_lu_to_s = Vec3::dist(lumen_closest_solid_vec,solid_vec);
    double dist_ti_to_s = Vec3::dist(tissue_closest_solid_vec,solid_vec);

    // ===== Get the surface normal from the closest point =====
    Vector_3 lumen_nor_vec;
    Vector_3 tissue_nor_vec;
    
    double l_n_pt[3],t_n_pt[3];
    
    for (int jj = 0; jj < 3; ++jj)
    {
      l_n_pt[jj] = *(lumen_wall_normal->GetTuple3(id_lumen_closest_solid_pt) + jj);
      t_n_pt[jj] = *(tissue_wall_normal->GetTuple3(id_tissue_closest_solid_pt) + jj);
    }

    lumen_nor_vec = l_n_pt;
    tissue_nor_vec = t_n_pt;

    lumen_nor_vec.normalize();
    tissue_nor_vec.normalize();

    // ===== Rescale the normals according to the weights =====
    lumen_nor_vec*=dist_lu_to_s/(dist_lu_to_s+dist_ti_to_s);
    tissue_nor_vec*=dist_ti_to_s/(dist_lu_to_s+dist_ti_to_s);

    // ===== Define radial normals by adding up ===== 
    Vector_3 loc_r_vec;

    if (lumen_nor_vec.dot_product(tissue_nor_vec) < 0)
      loc_r_vec = lumen_nor_vec - tissue_nor_vec;
    else
      loc_r_vec = lumen_nor_vec + tissue_nor_vec;
    loc_r_vec.normalize();

    loc_r_vec_array->InsertTuple3(ii, loc_r_vec.x(), loc_r_vec.y(), loc_r_vec.z());

    // ===== Define longitudinal normals ===== 
    Vector_3 temp_l_vec;

    for (size_t jj = 0; jj < Locators.size(); ++jj) 
    {
      // Get the branch centerline vtp from containers
      vtkSmartPointer<vtkPolyData> branch_poly = PolyDatas[jj];

      // Ensure the vtp is valid
      if (!branch_poly) {
          SYS_T::print_fatal_if( !branch_poly, "Error: SV_T:: Invalid branch centerline vtp from containers. \n");
          continue;
      }

      // Access the "radius_" data array
      vtkDataArray* dataArray = nullptr;
      for (int kk = 0; kk < branch_poly->GetPointData()->GetNumberOfArrays(); ++kk) {
          std::string arrayName = branch_poly->GetPointData()->GetArrayName(kk);
          if (arrayName.find("radius_") == 0) {
              dataArray = branch_poly->GetPointData()->GetArray(arrayName.c_str());
              break;
          }
      }

      if (!dataArray) {
          SYS_T::print_fatal_if( !dataArray, "Error: SV_T:: 'radius_' data array not found in dataset. \n");
          continue;
      }

      // ===== Define and get the closest point from solid point to the centerline point =====
      int id_centerline_closest_solid_pt=Locators[jj]->FindClosestPoint(solid_pt);
      double centerline_closest_solid_pt[3]; 
      branch_poly->GetPoint(id_centerline_closest_solid_pt, centerline_closest_solid_pt);

      // ===== Use Vector_3 store the closest point =====
      Vector_3 centerline_closest_solid_vec(centerline_closest_solid_pt[0],
                                            centerline_closest_solid_pt[1],
                                            centerline_closest_solid_pt[2]);

      // ===== Create splines with vtkParametricSpline ===== 
      vtkNew<vtkParametricSpline> spline;
      spline->SetPoints(branch_poly->GetPoints());

      // ===== Generate a geometric representation of the spline using vtkParametricFunctionSource ===== 
      vtkNew<vtkParametricFunctionSource> functionSource;
      functionSource->SetParametricFunction(spline);
      functionSource->Update();


      // ===== Compute the distance from solid point to the closest centerline =====
      double dist_center_to_s = Vec3::dist(centerline_closest_solid_vec,solid_vec);
      
      // ===== Update when the relative distance has a smaller value =====
      if ((dist_center_to_s / dataArray->GetComponent(id_centerline_closest_solid_pt, 0)) < dist_relevant[ii])
      {
        dist_relevant[ii] = dist_center_to_s / dataArray->GetComponent(id_centerline_closest_solid_pt, 0); 

        // ===== Computes the tangent vector at the given parameters ===== 
        double u[3], du[9];
        // ===== Parameterized value, ranging from 0 to 1, here we take the middle point ===== 
        u[0] = 0.5;  
        u[1] = u[2] = 0.0;

        spline->Evaluate(u, centerline_closest_solid_pt, du);
        Vector_3 temp(du[0],du[1],du[2]); 
        temp_l_vec = temp;
      }
    }

    Vector_3 loc_l_vec;
    Vector_3 loc_c_vec;

    // ===== Computes the circumferential direction vecters ===== 
    loc_c_vec = Vec3::cross_product(temp_l_vec, loc_r_vec);
    loc_c_vec.normalize();

    // ===== Computes the correct longitudinal direction vecters ===== 
    loc_l_vec = Vec3::cross_product(loc_r_vec, loc_c_vec);
    loc_l_vec.normalize();

    loc_l_vec_array->InsertTuple3(ii, loc_l_vec.x(), loc_l_vec.y(), loc_l_vec.z());
    loc_c_vec_array->InsertTuple3(ii, loc_c_vec.x(), loc_c_vec.y(), loc_c_vec.z());
    if (ii % 1000 == 0)
      printf("normal ii = %d\n", ii);
  }
  
  solid_grid->GetPointData()->AddArray(loc_r_vec_array);
  solid_grid->GetPointData()->AddArray(loc_l_vec_array);
  solid_grid->GetPointData()->AddArray(loc_c_vec_array);

  vtkNew<vtkXMLUnstructuredGridWriter> writer;
  std::string prefix_name_to_write = "out";
  prefix_name_to_write.append(solid_vtu_s.c_str());
  writer -> SetFileName( prefix_name_to_write.c_str() );
  writer -> SetInputData( solid_grid );
  writer -> Write();

  return 0;
}

// EOF
