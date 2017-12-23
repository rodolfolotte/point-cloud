#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_map.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Timer.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/Point_set_processing_3.h>

// IO
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/IO/write_ply_points.h>

// Simplification and smothing
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/estimate_scale.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>

// Outliers removal
#include <CGAL/remove_outliers.h>

// Normal and its orient
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/vcm_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>

// Structuring
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/regularize_planes.h>
#include <CGAL/structure_point_set.h>

#include <boost/foreach.hpp>

#include <utility>
#include <vector>
#include <string>
//#include <time.h>
//#include <cstdlib>
#include <iostream>
#include <ctime>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_set_3<Point> Points;
typedef std::pair<Point, Vector> PointAndNormal;
typedef std::vector<Point> PointList;
typedef std::vector<PointAndNormal> PointAndNormalList;
typedef CGAL::cpp11::array<double, 6> Covariance;
typedef CGAL::First_of_pair_property_map<PointAndNormal> Point_map;
typedef CGAL::Second_of_pair_property_map<PointAndNormal> Normal_map;

// Efficient RANSAC types
typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits<Kernel, PointAndNormalList, Point_map, Normal_map> Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits> Efficient_ransac;
typedef CGAL::Shape_detection_3::Plane<Traits> Plane;

#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

using namespace std;

std::string now()
{
  std::time_t t = std::time(0);
  char cstr[128];
  std::strftime(cstr, sizeof(cstr), "%d-%m-%Y %H:%M:%S", std::localtime(&t));

  return cstr;
}

std::string log(int option)
{
  string str;
  if (option == 0)
  {
    str = "[INFO | " + now() + "]";
  }
  else if (option == 1)
  {
    str = "[ERROR | " + now() + "]";
  }
  return str;
}

/* Read .OFF, .XYZ, .PWD, and PLY point cloud format files 
    input_filename: name of the file to be read */
Points readPCReturnsPoint_set_3(string input_filename)
{
  Points points;
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Open " << input_filename << " for reading..." << std::endl;
  string extension = input_filename.substr(input_filename.find_last_of('.'));
  std::ifstream stream(input_filename.c_str());
  bool success = false;

  if (!stream)
  {
    std::cerr << log(1) << " > Can't read input file " << std::endl;
  }

  if (extension == ".off" || extension == ".OFF")
  {
    success = CGAL::read_off_point_set(stream, points);
  }
  else if (extension == ".xyz" || extension == ".XYZ" || extension == ".txt" || extension == ".pwn" || extension == ".PWN")
  {
    success = CGAL::read_xyz_point_set(stream, points);
  }
  else if (extension == ".ply" || extension == ".PLY")
  {
    success = CGAL::read_ply_point_set(stream, points);
  }

  if (success)
  {
    int nb_points = points.size();
    std::cout << log(0) << " > File reading finished! Total: " << nb_points << " points, " << task_timer.time() << " seconds" << std::endl;
    task_timer.reset();

    if (nb_points == 0)
    {
      std::cerr << log(1) << " > Error: empty file" << std::endl;
    }
  }
  else
  {
    std::cerr << log(1) << " > Error: problems during the file reading" << std::endl;
  }

  return points;
}

/* Read .OFF, .XYZ, .PWD, and PLY point cloud format files  
    input_filename: name of the file to be read */
PointList readPCReturnsPointList(string input_filename)
{
  PointList points;
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Open " << input_filename << " for reading..." << std::endl;
  string extension = input_filename.substr(input_filename.find_last_of('.'));
  std::ifstream stream(input_filename.c_str());
  bool success = false;

  if (!stream)
  {
    std::cerr << log(1) << " > Can't read input file " << std::endl;
  }

  if (extension == ".off" || extension == ".OFF")
  {
    success = CGAL::read_off_points(stream, std::back_inserter(points), CGAL::Identity_property_map<Point>());
  }
  else if (extension == ".xyz" || extension == ".XYZ" || extension == ".txt" || extension == ".pwn" || extension == ".PWN")
  {
    success = CGAL::read_xyz_points(stream, std::back_inserter(points), CGAL::Identity_property_map<Point>());
  }
  else if (extension == ".ply" || extension == ".PLY")
  {
    success = CGAL::read_ply_points(stream, std::back_inserter(points), CGAL::Identity_property_map<Point>());
  }

  if (success)
  {
    int nb_points = points.size();
    std::cout << log(0) << " > File reading finished! Total: " << nb_points << " points, " << task_timer.time() << " seconds" << std::endl;
    task_timer.reset();

    if (nb_points == 0)
    {
      std::cerr << log(1) << " > Error: empty file" << std::endl;
    }
  }
  else
  {
    std::cerr << log(1) << " > Error: problems during the file reading" << std::endl;
  }

  return points;
}

/* Read .OFF, .XYZ, .PWD, and PLY point cloud format files  
    input_filename: name of the file to be read */
PointAndNormalList readPCReturnsPairList(string input_filename)
{
  PointAndNormalList points;
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Open " << input_filename << " for reading..." << std::endl;
  string extension = input_filename.substr(input_filename.find_last_of('.'));
  std::ifstream stream(input_filename.c_str());
  bool success = false;

  if (!stream)
  {
    std::cerr << log(1) << " > Can't read input file " << std::endl;
  }

  /*TODO: não seria read_XXX_points_and_normals ? */
  if (extension == ".off" || extension == ".OFF")
  {
    success = stream && CGAL::read_off_points_and_normals(stream, std::back_inserter(points), CGAL::First_of_pair_property_map<PointAndNormal>(), CGAL::Second_of_pair_property_map<PointAndNormal>());
  }
  else if (extension == ".xyz" || extension == ".XYZ" || extension == ".txt" || extension == ".pwn" || extension == ".PWN")
  {
    success = CGAL::read_xyz_points_and_normals(stream, std::back_inserter(points), CGAL::First_of_pair_property_map<PointAndNormal>(), CGAL::Second_of_pair_property_map<PointAndNormal>());
  }
  else if (extension == ".ply" || extension == ".PLY")
  {
    success = CGAL::read_ply_points_and_normals(stream, std::back_inserter(points), CGAL::First_of_pair_property_map<PointAndNormal>(), CGAL::Second_of_pair_property_map<PointAndNormal>());
  }

  if (success)
  {
    int nb_points = points.size();
    std::cout << log(0) << " > File reading finished! Total: " << nb_points << " points, " << task_timer.time() << " seconds" << std::endl;
    task_timer.reset();

    if (nb_points == 0)
    {
      std::cerr << log(1) << " > Error: empty file" << std::endl;
    }
  }
  else
  {
    std::cerr << log(1) << " > Error: problems during the file reading" << std::endl;
  }

  return points;
}

/* Write .OFF, .XYZ, .PWD, and PLY point cloud format files  
    points: point cloud Point_set_3 typed to be save
    output_filename: name of the file to be save */
void writePoint_set_3(Points points, string output_filename)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Writing file " << output_filename << " with " << points.size() << " points... " << std::endl;
  std::string extension = output_filename.substr(output_filename.find_last_of('.'));
  ofstream stream(output_filename.c_str());
  bool success = false;

  if (extension == ".xyz" || extension == ".XYZ" || extension == ".pwn" || extension == ".PWN")
  {
    success = CGAL::write_xyz_point_set(stream, points);
  }
  else if (extension == ".off" || extension == ".OFF")
  {
    success = CGAL::write_off_point_set(stream, points);
  }
  else if (extension == ".ply" || extension == ".PLY")
  {
    success = CGAL::write_ply_point_set(stream, points);
  }

  if (success)
  {
    std::cout << log(0) << " > " << output_filename << " successfull saved," << task_timer.time() << " seconds" << std::endl;
  }
  else
  {
    std::cerr << log(1) << " > Error: cannot write file " << output_filename << std::endl;
  }
}

/* Write .OFF, .XYZ, .PWD, and PLY point cloud format files   
    points: point cloud std::vector<Point> typed to be save
    output_filename: name of the file to be save */
void writePointList(PointList points, string output_filename)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Writing file " << output_filename << " with " << points.size() << " points... " << std::endl;
  std::string extension = output_filename.substr(output_filename.find_last_of('.'));
  ofstream stream(output_filename.c_str());
  bool success = false;

  if (extension == ".xyz" || extension == ".XYZ" || extension == ".pwn" || extension == ".PWN")
  {
    success = CGAL::write_xyz_points(stream, points.begin(), points.end());
  }
  else if (extension == ".off" || extension == ".OFF")
  {
    success = CGAL::write_off_points(stream, points.begin(), points.end());
  }
  else if (extension == ".ply" || extension == ".PLY")
  {
    success = CGAL::write_ply_points(stream, points.begin(), points.end());
  }

  if (success)
  {
    std::cout << log(0) << " > " << output_filename << " successfull saved, " << task_timer.time() << " seconds" << std::endl;
  }
  else
  {
    std::cerr << log(1) << " > Error: cannot write file " << output_filename << std::endl;
  }
}

/* Write .OFF, .XYZ, .PWD, and PLY point cloud format files   
    points: point cloud std::vector<Point> typed to be save
    output_filename: name of the file to be save */
void writePairList(PointAndNormalList points, string output_filename)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Writing file " << output_filename << " with " << points.size() << " points... " << std::endl;
  std::string extension = output_filename.substr(output_filename.find_last_of('.'));
  ofstream stream(output_filename.c_str());
  bool success = false;

  if (extension == ".xyz" || extension == ".XYZ" || extension == ".pwn" || extension == ".PWN")
  {
    success = CGAL::write_xyz_points_and_normals(stream, points.begin(), points.end(),
                                                 CGAL::First_of_pair_property_map<PointAndNormal>(),
                                                 CGAL::Second_of_pair_property_map<PointAndNormal>());
  }
  else if (extension == ".off" || extension == ".OFF")
  {
    success = CGAL::write_off_points_and_normals(stream, points.begin(), points.end(),
                                                 CGAL::First_of_pair_property_map<PointAndNormal>(),
                                                 CGAL::Second_of_pair_property_map<PointAndNormal>());
  }
  else if (extension == ".ply" || extension == ".PLY")
  {
    success = CGAL::write_ply_points_and_normals(stream, points.begin(), points.end(),
                                                 CGAL::First_of_pair_property_map<PointAndNormal>(),
                                                 CGAL::Second_of_pair_property_map<PointAndNormal>());
  }

  if (success)
  {
    std::cout << log(0) << " > " << output_filename << " successfull saved, " << task_timer.time() << " seconds" << std::endl;
  }
  else
  {
    std::cerr << log(1) << " > Error: cannot write file " << output_filename << std::endl;
  }
}

/* Simplify point set 
    points: point cloud (vector<Point_3>) */
void simplify(PointList &points)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Grid simplifying..." << std::endl;
  std::cout << log(0) << " > Estimating global range scale..." << std::endl;
  FT range_scale = CGAL::estimate_global_range_scale(points.begin(), points.end());

  std::cout << log(0) << " > Simplifying grid with range scale of " << range_scale << "..." << std::endl;
  points.erase(CGAL::grid_simplify_point_set(points.begin(), points.end(), range_scale), points.end());

  std::cout << log(0) << " > Point cloud simplified. Total: " << points.size() << " points, " << task_timer.time() << " seconds" << std::endl;
}

void simplifyWLOP(PointList &points, double retain_percentage, double neighbor_radius)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Grid simplifying with WLOP algorithm: Percentage to retain " << retain_percentage << "%, and neighbor radius " << neighbor_radius << " ..." << std::endl;

  PointList output;
  CGAL::wlop_simplify_and_regularize_point_set<Concurrency_tag>(points.begin(), points.end(), std::back_inserter(output), retain_percentage, neighbor_radius);
  points = output;

  std::cout << log(0) << " > Point cloud simplified with WLOP algorithm. Total: " << points.size() << " points, " << task_timer.time() << " seconds" << std::endl;
}

/* Simplify point set 
    points: point cloud (vector<Point_3>) */
void simplify(PointList &points, double gridSize)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Grid simplifying..." << std::endl;
  points.erase(CGAL::grid_simplify_point_set(points.begin(), points.end(), gridSize), points.end());
  std::cout << log(0) << " > Point cloud simplified. Total: " << points.size() << " points, " << task_timer.time() << " seconds" << std::endl;
}

/* Jet Smothing point cloud - use estimated k as scale for jet smoothing and range for grid simplification
    points: point cloud */
void smooth(PointList &points)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Smoothing the point cloud..." << std::endl;
  std::cout << log(0) << " > Estimating global k neighbor scale..." << std::endl;
  std::size_t k_scale = CGAL::estimate_global_k_neighbor_scale(points.begin(), points.end());

  std::cout << log(0) << " > Smothing point cloud by Jet algorithm with k scale of " << k_scale << "..." << std::endl;
  CGAL::jet_smooth_point_set<Concurrency_tag>(points.begin(), points.end(), static_cast<unsigned int>(k_scale));
  std::cout << log(0) << " > Point cloud smoothed. Total: " << points.size() << " points, " << task_timer.time() << " seconds" << std::endl;
}

/* Computes normals direction by Principal Component Analysis - input points + output normals */
void estimateNormalsWithPCA(PointAndNormalList &points, unsigned int nb_neighbors)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Estimates Normals Direction by PCA with neighbors = " << nb_neighbors << "..." << std::endl;

  // Estimates normals direction.
  // Note: pca_estimate_normals() requires an iterator over points
  // as well as property maps to access each point's position and normal.
  CGAL::pca_estimate_normals<Concurrency_tag>(points.begin(), points.end(),
                                              CGAL::First_of_pair_property_map<PointAndNormal>(),
                                              CGAL::Second_of_pair_property_map<PointAndNormal>(),
                                              nb_neighbors);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cout << log(0) << " > Estimate Normals Direction by PCA concluded, " << points.size() << " points, and " << (memory >> 20) << " Mb allocated, " << task_timer.time() << " seconds" << std::endl;
}

/* Computes normals direction by Jet Fitting - input points + output normals
    points: point cloud
    neighbors: threshold of neighbors */
void estimateNormalsWithJET(Points &points, int neighbors)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Estimates Normals Direction by Jet Fitting (neighbors = " << neighbors << ")..." << std::endl;
  CGAL::jet_estimate_normals<CGAL::Sequential_tag>(points, neighbors);
  long memory = CGAL::Memory_sizer().virtual_size();
  std::cout << log(0) << " > Estimate Normals Direction by Jet Fitting concluded, " << (memory >> 20) << " Mb allocated, " << task_timer.time() << " seconds" << std::endl;
}

/* Computes normals direction by Jet Fitting - input points + output normals */
void estimateNormalsWithJET(PointAndNormalList &points, unsigned int neighbors)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Estimates Normals Direction by Jet Fitting (k=" << neighbors << ")..." << std::endl;

  // Estimates normals direction.
  // Note: jet_estimate_normals() requires an iterator over points
  // + property maps to access each point's position and normal.
  CGAL::jet_estimate_normals<Concurrency_tag>(points.begin(), points.end(),
                                              CGAL::First_of_pair_property_map<PointAndNormal>(),
                                              CGAL::Second_of_pair_property_map<PointAndNormal>(),
                                              neighbors);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cout << log(0) << " > Estimate Normals Direction by Jet Fitting concluded, " << (memory >> 20) << " Mb allocated, " << task_timer.time() << " seconds" << std::endl;
}

/* Compute normals direction using the VCM
   input: point cloud
   R: radius of the offset
   r: radius used during the convolution */
void estimateNormalsWithVCM(PointAndNormalList &points, double R, double r)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Estimates Normals Direction using VCM (R=" << R << " and r=" << r << ")..." << std::endl;

  // Estimates normals direction.
  // Note: vcm_estimate_normals() requires an iterator over points
  // + property maps to access each point's position and normal.
  CGAL::vcm_estimate_normals(points.begin(), points.end(),
                             CGAL::First_of_pair_property_map<PointAndNormal>(),
                             CGAL::Second_of_pair_property_map<PointAndNormal>(),
                             R, r);

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cout << log(0) << " > Estimate Normals Direction by VCM concluded, " << (memory >> 20) << " Mb allocated, " << task_timer.time() << " seconds" << std::endl;
}

/* Normal orientation using a Minimum Spanning Tree. - input points + input/output normals 
   @book{hoppe1992,
    title={Surface reconstruction from unorganized points},
    author={Hoppe, Hugues and DeRose, Tony and Duchamp, Tom and McDonald, John and Stuetzle, Werner},
    volume={26},
    number={2},
    year={1992},
    publisher={ACM}
  } */
void orientNormalsWithMST(PointAndNormalList &points, unsigned int nb_neighbors_mst)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Orients Normals with a Minimum Spanning Tree (k=" << nb_neighbors_mst << ")..." << std::endl;

  // Orients normals.
  // Note: mst_orient_normals() requires an iterator over points
  // as well as property maps to access each point's position and normal.
  PointAndNormalList::iterator unoriented_points_begin =
      CGAL::mst_orient_normals(points.begin(), points.end(),
                               CGAL::First_of_pair_property_map<PointAndNormal>(),
                               CGAL::Second_of_pair_property_map<PointAndNormal>(),
                               nb_neighbors_mst);

  // Optional: delete points with an unoriented normal
  // if you plan to call a reconstruction algorithm that expects oriented normals.
  points.erase(unoriented_points_begin, points.end());

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cout << log(0) << " > Orient Normals by MST concluded, " << (memory >> 20) << " Mb allocated, " << task_timer.time() << " seconds" << std::endl;
}

/* Removes outliers using erase-remove idiom
    option: 0: I KNOW the ratio of outliers present in the point set
            1: I DON'T KNOW the ratio of outliers present in the point set
    nb_neighbors: considers nearest neighbor points
    limitToRemove: 100.0 = No limit on the number of outliers to remove
    spacingAverage: Point with distance above 2*average_spacing are considered outliers */
void removeOutliers(PointList &points, int option, int nb_neighbors, double limitToRemove, double spacingAvarege)
{
  CGAL::Timer task_timer;
  task_timer.start();

  // Estimate scale of the point set with average spacing
  std::cout << log(0) << " Removing outliers..." << std::endl;
  std::cout << log(0) << " > Computing average spacing..." << std::endl;
  const double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(points.begin(), points.end(), nb_neighbors);
  std::cout << log(0) << " > Removing outliers with " << average_spacing << " of average spacing..." << std::endl;

  if (option == 0)
  {
    std::vector<Point>::iterator first_to_remove = CGAL::remove_outliers(points.begin(), points.end(), CGAL::Identity_property_map<Point>(), nb_neighbors,
                                                                         limitToRemove, spacingAvarege * average_spacing);

    std::cout << log(0) << " > " << (limitToRemove * std::distance(first_to_remove, points.end()) / (double)(points.size()))
              << "% of the points are considered outliers when using a distance threshold of "
              << spacingAvarege * average_spacing << ", " << task_timer.time() << " seconds" << std::endl;
  }
  else
  {
    points.erase(CGAL::remove_outliers(points.begin(), points.end(),
                                       CGAL::Identity_property_map<Point>(),
                                       nb_neighbors,
                                       limitToRemove,
                                       0.0),
                 points.end());

    // Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
    std::vector<Point>(points).swap(points);
  }

  std::cout << log(0) << " > Outliers removal concluded. Total: " << points.size() << " points, " << task_timer.time() << " seconds" << std::endl;
}

/* TODO */
void structurePointSet(PointAndNormalList &points, double epsilon, bool parallelism, bool orthogonality, bool coplanarity, bool zSymmetry, int degreeOfTolerance)
{
  CGAL::Timer task_timer;
  task_timer.start();

  std::cout << log(0) << " Structuring the point cloud..." << std::endl;
  std::cout << log(0) << " > Detecting planes with RANSAC..." << std::endl;
  Efficient_ransac ransac;
  ransac.set_input(points);
  ransac.add_shape_factory<Plane>();
  ransac.detect();

  /* regularize parallelism? 
     regularize orthogonality? 
     regularize coplanarity? 
     regularize Z-symmetry? 
     degrees of tolerance for parallelism/orthogonality */
  std::cout << log(0) << " > Regularizing planes with parallelism[" << parallelism << "], orthogonality[" << orthogonality << "], coplanarity[" << coplanarity << "], zSymmetry[" << zSymmetry << "], and degree of tolerance: " << degreeOfTolerance << "..." << std::endl;
  CGAL::regularize_planes(ransac, parallelism, orthogonality, coplanarity, zSymmetry, degreeOfTolerance);

  std::cout << log(0) << " > Structuring the point set with epsilon " << epsilon << " ..." << std::endl;
  CGAL::structure_point_set(points.begin(), points.end(), std::back_inserter(points), ransac, epsilon); // epsilon for structuring points

  std::cout << log(0) << " > Structured point set generated, " << points.size() << " points, " << task_timer.time() << " seconds" << std::endl;
}

bool str2bool(std::string str)
{
  if (str == "true")
  {
    return true;
  }
  else
  {
    return false;
  }
}

int main(int argc, char *argv[])
{
  string input_filename = argv[1];
  string output_filename = argv[2];

  string path = output_filename.substr(0, output_filename.find_last_of("\\/"));
  string log_filename = path + "/_log.txt";
  FILE *fileOut = freopen(log_filename.c_str(), "a", stdout);
  FILE *fileErr = freopen(log_filename.c_str(), "a", stderr);

  std::string procedure = "-simplify";
  std::string param_procedure = "wlop";

  unsigned int nb_neighbors_pca_normals = 18;
  unsigned int nb_neighbors_jet_fitting_normals = 18;
  unsigned int nb_neighbors_mst = 18;

  double offset_radius_vcm = 0.2;
  double convolve_radius_vcm = 0.1;

  double neighbor_radius = 0.5;
  double retain_percentage = 2.0;

  double epsilon = 0.015;
  bool parallelism = true;
  bool orthogonality = true;
  bool coplanarity = false;
  bool zSymmetry = true;
  int degreeOfTolerance = 10;

  CGAL::Timer task_timer;
  Points points;
  PointList pointList;
  PointAndNormalList pairs;

  std::cout << " ---- 3D Reconstruction: A facade-feature extractor " << std::endl;

  if (argc - 1 < 2)
  {
    std::cout << " Prepocessing point sets." << std::endl;
    std::cout << " " << std::endl;
    std::cout << " Usage: " << argv[0] << " input output [procedure] [method] [n_param] [n_param_value]" << std::endl;
    std::cout << " Example:          " << argv[0] << " input output -simplify random " << std::endl;
    std::cout << " (default values): " << argv[0] << " input output -simplify wlop " << std::endl;
    std::cout << " (custom values)   " << argv[0] << " input output -simplify wlop -neighbor_radius 0.8 -retain_percentage 0.5" << std::endl;
    std::cout << "                   " << argv[0] << " input output -normal plane -nb_neighbors_pca 15" << std::endl;
    std::cout << "                   " << argv[0] << " input output -orient -nb_neighbors_mst 12" << std::endl;
    std::cout << " Input file formats are .off, .xyz .pwn and .ply" << std::endl;
    std::cout << " Output file formats are .off, .xyz, and .ply." << std::endl;
    std::cout << " Options:" << std::endl;
    std::cout << "  -simplify wlop|random|hierarchy|var_max               Simplify the point set based on its neighbor and percent retaintion (default=wlop)" << std::endl;
    std::cout << "    -nb_wlop_simplify <double>                           Number of neighbors to compute tangent plane (default=0.5)" << std::endl;
    std::cout << "    -retain_percentage <double>                         Percentage of points to be reduced (default=0.2)" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "  -smooth                                               Smooth point set based on ... REF ... algorithm" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "  -outlier                                              Remove outliers based on ... REF ... algorithm" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "  -normal plane|quadric|vcm                             Estimates normals direction and orientation using a tangent plane or quadric or vcm (default=plane)" << std::endl;
    std::cout << "    -nb_neighbors_pca <int>                             Number of neighbors to compute tangent plane (default=18)" << std::endl;
    std::cout << "    -nb_neighbors_jet_fitting <int>                     Number of neighbors to compute quadric (default=18)" << std::endl;
    std::cout << "    -offset_radius_vcm <double>                         Offset radius to compute VCM (default=0.1)" << std::endl;
    std::cout << "    -convolve_radius_vcm <double>                       Convolve radius to compute VCM (default=0.1)" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "  -orient                                               Orient normals using a Minimum Spanning Tree" << std::endl;
    std::cout << "    -nb_neighbors_mst <int>                             Number of neighbors to compute the MST (default=18)" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "  -structure                                            Surface reconstruction through point set structuring" << std::endl;
    std::cout << "    -epsilon <double>                                   Epsilon ... (default=0.015)" << std::endl;
    std::cout << "    -paral <bool>                                       Regularize parallelism ? (default=true)" << std::endl;
    std::cout << "    -ortho <bool>                                       Regularize orthogonality ? (default=true)" << std::endl;
    std::cout << "    -coplan <bool>                                      Regularize coplanarity ? (default=false)" << std::endl;
    std::cout << "    -zsymmetry <bool>                                   Regularize z-symmetry ? (default=true)" << std::endl;
    std::cout << "    -degree <int>                                       Degree of tolerance (default=10)" << std::endl;
    std::cout << " " << std::endl;

    return EXIT_FAILURE;
  }

  task_timer.start();

  for (int i = 3; i + 1 < argc; ++i)
  {
    procedure = argv[i];
    if (procedure == "-simplify")
    {
      param_procedure = argv[++i];
      pointList = readPCReturnsPointList(input_filename);

      if (param_procedure == "wlop")
      {        
        if (std::string(argv[++i]) == "-nb_wlop_simplify")
          neighbor_radius = atof(argv[++i]);
        
        if (std::string(argv[++i]) == "-retain_percentage")
          retain_percentage = atof(argv[++i]);
        
        simplifyWLOP(pointList, retain_percentage, neighbor_radius);
        //simplify(pointList, simplification_factor);
        
        //writePoint_set_3(points, output_filename);
        writePointList(pointList, output_filename);
      }
      else if (param_procedure == "random")
      {
        simplify(pointList);
        writePointList(pointList, output_filename);
      }
      else if (param_procedure == "hierarchy")
      {
        //TODO
      }
      else if (param_procedure == "var_max")
      {
        //TODO
      }
      else
      {
        std::cerr << log(1) << " > Invalid parameter for simplifing point set! Check it and try again!" << std::endl;
      }
    }
    else if (procedure == "-smooth")
    {
      //smooth(points);
      //smooth(pointList);
      //TODO
      std::cerr << log(1) << " > Sorry! Nothing has been done to smooth your point cloud!" << argv[i] << "" << std::endl;
    }
    else if (procedure == "-outlier")
    {
      //TODO
      //removeOutliers(points, 0, 12, 5.0, 2.0);
      //removeOutliers(pointList, 0, neighbor_factor, limitToRemove_factor, spacingAverage_factor);
      std::cerr << log(1) << " > Sorry! Nothing has been done to remove outlier from your point cloud!" << argv[i] << "" << std::endl;
    }
    else if (procedure == "-normal")
    {
      param_procedure = argv[++i];
      pairs = readPCReturnsPairList(input_filename);

      if (param_procedure == "plane")
      {
        if (std::string(argv[++i]) == "-nb_neighbors_pca")
          nb_neighbors_pca_normals = atoi(argv[++i]);

        estimateNormalsWithPCA(pairs, nb_neighbors_pca_normals);
        writePairList(pairs, output_filename);
      }
      else if (param_procedure == "quadric")
      {
        if (std::string(argv[++i]) == "-nb_neighbors_jet_fitting")
          nb_neighbors_jet_fitting_normals = atoi(argv[++i]);

        estimateNormalsWithJET(pairs, nb_neighbors_jet_fitting_normals);
        writePairList(pairs, output_filename);
      }
      else if (param_procedure == "vcm")
      {

        if (std::string(argv[++i]) == "-offset_radius_vcm")
          offset_radius_vcm = atof(argv[++i]);

        if (std::string(argv[++i]) == "-convolve_radius_vcm")
          convolve_radius_vcm = atof(argv[++i]);

        estimateNormalsWithVCM(pairs, offset_radius_vcm, convolve_radius_vcm);
        writePairList(pairs, output_filename);
      }
      else
      {
        std::cerr << log(1) << " > Invalid parameter for estimating normals! Check it and try again!" << std::endl;
      }
    }
    else if (procedure == "-orient")
    {
      pairs = readPCReturnsPairList(input_filename);

      if (std::string(argv[++i]) == "-nb_neighbors_mst")
      {
        nb_neighbors_mst = atoi(argv[++i]);

        orientNormalsWithMST(pairs, nb_neighbors_mst);
        writePairList(pairs, output_filename);
      }
      else
      {
        std::cerr << log(1) << " > Invalid parameter for orienting normal! Check it and try again!" << std::endl;
      }
    }
    else if (procedure == "-structure")
    {
      pairs = readPCReturnsPairList(input_filename);

      if (std::string(argv[++i]) == "-epsilon")
        epsilon = atof(argv[++i]);

      if (std::string(argv[++i]) == "-paral")
        parallelism = str2bool(std::string(argv[++i]));

      if (std::string(argv[++i]) == "-ortho")
        orthogonality = str2bool(std::string(argv[++i]));

      if (std::string(argv[++i]) == "-coplan")
        coplanarity = str2bool(std::string(argv[++i]));

      if (std::string(argv[++i]) == "-zsymmetry")
        zSymmetry = str2bool(std::string(argv[++i]));

      if (std::string(argv[++i]) == "-degree")
        degreeOfTolerance = atoi(argv[++i]);

      structurePointSet(pairs, epsilon, parallelism, orthogonality, coplanarity, zSymmetry, degreeOfTolerance);
      writePairList(pairs, output_filename);
    }
    else
    {
      std::cerr << log(1) << " > Invalid option " << argv[i] << "" << std::endl;
    }
  }

  std::cout << log(0) << " End of process! Total: " << task_timer.time() << " seconds" << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;
}