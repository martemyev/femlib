#include "parameters.h"
#include "config.h"
#include "auxiliary_functions.h"
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include <iostream>

namespace po = boost::program_options;



Parameters::Parameters(int argc, char **argv)
{
  default_parameters(); // initialize all parameters
  if (argc > 1)
    read_from_command_line(argc, argv); // change some (or all) parameters from default to custom ones
}



void Parameters::default_parameters()
{
  MESH_DIR = HOME_DIRECTORY + "/projects/tat_gmsfem/brandnew/meshes/";
  GEO_DIR  = HOME_DIRECTORY + "/projects/tat_gmsfem/brandnew/geometries/";

  RES_TOP_DIR = "../results/"; // this top level directory containing all results usually exists on the same level as 'build', 'sources', 'headers' directories
  RES_DIR = ""; // should be changed and based on some parameters
  VTU_DIR = "vtu/"; // should be added to RES_DIR after generating of the latter
  SOL_DIR = "sol/"; // should be added to RES_DIR after generating of the latter
  TIME_FILE = "time.txt"; // should be added to RES_DIR after generating of the latter
  INFO_FILE = "info.txt"; // should be added to RES_DIR after generating of the latter

  MESH_FILE = "mesh.msh";  // should be added to MESH_DIR after establishing of the latter (means that MESH_DIR can be changed from parameter file of command line)

  TIME_SCHEME = EXPLICIT;
  X_BEG = Y_BEG = 0.;
  X_END = Y_END = 1.;
  //CL = 0.1;
  TIME_BEG = 0.;
  TIME_END = 1.;
  TIME_STEP = 1.;
  N_TIME_STEPS = 1;
  FE_ORDER = 1;
  //QUAD_ORDER = 3;
  SOURCE_FREQUENCY = 20;
  SOURCE_SUPPORT = 10;
  SOURCE_CENTER_X = 0.5;
  SOURCE_CENTER_Y = 0.5;

  N_SUBDOMAINS = 1; // by default there is only one domain with homogeneous coefficients
  COEF_A_FILES.push_back("");
  COEF_B_FILES.push_back("");
  COEF_A_VALUES.push_back(1.);
  COEF_B_VALUES.push_back(1.);

  IPDG_GAMMA = 2.0;

  PRINT_VTU = 0; // don't print .vtu files by default
  SAVE_SOL = 0; // don't print .sol files by default
  PRINT_INFO = 0; // don't print an information to console on each time step
  VTU_STEP = 1; // print the .vtu file on each time step
  SOL_STEP = 1; // save the .dat file with solution on each time step

  N_COARSE_ELEMENTS = 1;
}



void Parameters::read_from_command_line(int argc, char **argv)
{
  std::string time_scheme = (TIME_SCHEME == EXPLICIT ? "explicit" : "crank-nicolson");
//  std::string coef_a_files = "";
//  std::string coef_b_files = "";
//  std::string coef_a_values = "";
//  std::string coef_b_values = "";

//  for (coef_a

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("meshfile", po::value<std::string>(),  std::string("name of mesh file (" + MESH_FILE + ")").c_str())
    ("meshdir",  po::value<std::string>(),  std::string("path to a directory with meshes (" + MESH_DIR + ")").c_str())
    ("scheme",   po::value<std::string>(),  std::string("time scheme (" + time_scheme + ")").c_str())
    ("tend",     po::value<double>(),       std::string("time ending (" + d2s(TIME_END) + ")").c_str())
    ("tstep",    po::value<double>(),       std::string("time step (" + d2s(TIME_STEP) + ")").c_str())
    ("nt",       po::value<unsigned int>(), std::string("number of time steps (" + d2s(N_TIME_STEPS) + ")").c_str())
    ("fe",       po::value<unsigned int>(), std::string("order of fe basis functions (" + d2s(FE_ORDER) + ")").c_str())
    ("nsub",     po::value<unsigned int>(), std::string("number of subdomains with different physical properties (" + d2s(N_SUBDOMAINS) + ")").c_str())
    ("aXfile",   po::value<std::string>(),  std::string("name of the file with coefficient alpha distribution in domain number X [1,..,nsub] ()").c_str())
    ("bXfile",   po::value<std::string>(),  std::string("name of the file with coefficient beta distribution in domain number X [1,..,nsub] ()").c_str())
    ("aXval",    po::value<double>(),       std::string("the value of coefficient alpha in domain number X [1,..,nsub] (" + d2s(COEF_A_VALUES[0]) + ")").c_str())
    ("bXval",    po::value<double>(),       std::string("the value of coefficient beta in domain number X [1,..,nsub] (" + d2s(COEF_B_VALUES[0]) + ")").c_str())
    ("vtu",      po::value<bool>(),         std::string("whether we need to print .vtu files (" + d2s(PRINT_VTU) + ")").c_str())
    ("sol",      po::value<bool>(),         std::string("whether we need to save .dat files with solutions (" + d2s(SAVE_SOL) + ")").c_str())
    ("inf",      po::value<bool>(),         std::string("whether we need to print some info during calculations (" + d2s(PRINT_INFO) + ")").c_str())
    ("vtu_step", po::value<unsigned int>(), std::string("if we need to print .vtu files then how often. every (vtu_step)-th file will be printed (" + d2s(VTU_STEP) + ")").c_str())
    ("sol_step", po::value<unsigned int>(), std::string("if we need to save .dat files then how often. every (sol_step)-th file will be saved (" + d2s(SOL_STEP) + ")").c_str())
    ("x1",       po::value<double>(),       std::string("X_END (" + d2s(X_END) + ")").c_str())
    ("y1",       po::value<double>(),       std::string("Y_END (" + d2s(Y_END) + ")").c_str())
    ("gamma",    po::value<double>(),       std::string("parameter of the IPDG method (" + d2s(IPDG_GAMMA) + ")").c_str())
    ("f0",       po::value<double>(),       std::string("source frequency (" + d2s(SOURCE_FREQUENCY) + ")").c_str())
    ("p",        po::value<double>(),       std::string("source parameter (" + d2s(SOURCE_SUPPORT) + ")").c_str())
    ("xcen",     po::value<double>(),       std::string("x-center of the source (" + d2s(SOURCE_CENTER_X) + ")").c_str())
    ("ycen",     po::value<double>(),       std::string("y-center of the source (" + d2s(SOURCE_CENTER_Y) + ")").c_str())
  ;

  po::variables_map vm;
  po::parsed_options parsed = po::command_line_parser(argc, argv).options(desc).allow_unregistered().run();
  //po::store(po::parse_command_line(argc, argv, desc), vm);
  po::store(parsed, vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << "\n";
    exit(0);
  }

  if (vm.count("meshfile"))
    MESH_FILE = vm["meshfile"].as<std::string>();

  if (vm.count("meshdir"))
    MESH_DIR = vm["meshdir"].as<std::string>();

  if (vm.count("scheme"))
  {
    std::string scheme_name = vm["scheme"].as<std::string>();
    if (scheme_name == "explicit")
      TIME_SCHEME = EXPLICIT;
    else if (scheme_name == "cn" || scheme_name == "crank-nicolson")
      TIME_SCHEME = CRANK_NICOLSON;
    else
      require(false, "Unknown time scheme : " + scheme_name);
  }

  require(!(vm.count("tstep") && vm.count("nt") && vm.count("tend")),
          "tstep, nt and tend parameters cannot be used together - maximum two of them");

  if (vm.count("tend"))
    TIME_END = vm["tend"].as<double>();
  if (vm.count("tstep"))
    TIME_STEP = vm["tstep"].as<double>();
  if (vm.count("nt"))
    N_TIME_STEPS = vm["nt"].as<unsigned int>();

  if (!vm.count("nt"))
    N_TIME_STEPS = int((TIME_END - TIME_BEG) / TIME_STEP);
  if (!vm.count("tstep"))
    TIME_STEP = (TIME_END - TIME_BEG) / N_TIME_STEPS;
  if (!vm.count("tend"))
    TIME_END = TIME_BEG + N_TIME_STEPS * TIME_STEP;

  require(TIME_BEG >= 0., "time of begin cannot be negative");
  require(TIME_END > TIME_BEG, "time of end is not bigger than time of begin");
  require(TIME_STEP <= (TIME_END - TIME_BEG),
          "time step cannot be bigger than the whole time: time_step = " + d2s(TIME_STEP) + ", total_time = " + d2s(TIME_END - TIME_BEG));
  require(N_TIME_STEPS > 0, "the number of time steps cannot be less than 1");
  require(fabs((TIME_END - TIME_BEG - N_TIME_STEPS * TIME_STEP) / TIME_END) < 1e-14,
          "time parameters (TIME_END, TIME_STEP and N_TIME_STEPS) do not conform to each other");

  if (vm.count("fe"))
    FE_ORDER = vm["fe"].as<unsigned int>();
  require(FE_ORDER == 1, "This order of basis functions (" + d2s(FE_ORDER) + ") is not implemented");

  if (vm.count("vtu"))
    PRINT_VTU = vm["vtu"].as<bool>();
  if (vm.count("sol"))
    SAVE_SOL = vm["sol"].as<bool>();
  if (vm.count("inf"))
    PRINT_INFO = vm["inf"].as<bool>();
  if (vm.count("vtu_step"))
    VTU_STEP = vm["vtu_step"].as<unsigned int>();
  if (vm.count("sol_step"))
    SOL_STEP = vm["sol_step"].as<unsigned int>();

  if (vm.count("x1"))
    X_END = vm["x1"].as<double>();
  if (vm.count("y1"))
    Y_END = vm["y1"].as<double>();
  require(X_END > X_BEG && Y_END > Y_BEG, "Wrong limits of the computational domain");

  if (vm.count("gamma"))
    IPDG_GAMMA = vm["gamma"].as<double>();
  require(IPDG_GAMMA >= 0, "Gamma is wrong");

  if (vm.count("f0"))
    SOURCE_FREQUENCY = vm["f0"].as<double>();
  require(SOURCE_FREQUENCY > 0, "Frequency is wrong");
  if (vm.count("p"))
    SOURCE_SUPPORT = vm["p"].as<double>();
  require(SOURCE_SUPPORT > 0, "Source support is wrong");
  if (vm.count("xcen"))
    SOURCE_CENTER_X = vm["xcen"].as<double>();
  if (vm.count("ycen"))
    SOURCE_CENTER_Y = vm["ycen"].as<double>();

  if (vm.count("nsub"))
    N_SUBDOMAINS = vm["nsub"].as<unsigned int>();
  require(N_SUBDOMAINS < _N_SUBDOMAINS_LIMIT, "Too many subdomains");

  COEF_A_FILES.resize(N_SUBDOMAINS, ""); // allocate memory for file names for each subdomain
  COEF_B_FILES.resize(N_SUBDOMAINS, ""); // allocate memory for file names for each subdomain
  COEF_A_VALUES.resize(N_SUBDOMAINS, 1); // allocate memory for values for each subdomain
  COEF_B_VALUES.resize(N_SUBDOMAINS, 1); // allocate memory for values for each subdomain

  for (int i = 1; i <= N_SUBDOMAINS; ++i)
  {
    std::string param_name_afile = "a" + d2s(i) + "file";
    std::string param_name_bfile = "b" + d2s(i) + "file";
    std::string param_name_aval  = "a" + d2s(i) + "val";
    std::string param_name_bval  = "b" + d2s(i) + "val";

    unsigned int count_a = 0; // counter for values of coef alpha for each subdomain
    if (vm.count(param_name_afile.c_str()))
    {
      ++count_a;
      COEF_A_FILES[i] = vm[param_name_afile.c_str()].as<std::string>();
    }
    if (vm.count(param_name_aval.c_str()))
    {
      ++count_a;
      COEF_A_VALUES[i] = vm[param_name_aval.c_str()].as<double>();
    }
    require(count_a != 2, "error");
    // if there is no mention of coefficient alpha in this domain,
    // coefficient alpha in this domain is equal to 1 (default value)

    unsigned int count_b = 0; // counter for values of coef beta for each subdomain
    if (vm.count(param_name_bfile.c_str()))
    {
      ++count_b;
      COEF_B_FILES[i] = vm[param_name_bfile.c_str()].as<std::string>();
    }
    if (vm.count(param_name_bval.c_str()))
    {
      ++count_b;
      COEF_B_VALUES[i] = vm[param_name_bval.c_str()].as<double>();
    }
    require(count_b != 2, "error");
    // if there is no mention of coefficient beta in this domain,
    // coefficient beta in this domain is equal to 1 (default value)
  }

}



Parameters::~Parameters()
{ }



std::string Parameters::print() const
{
  std::string time_scheme_name[] = { "explicit",
                                     "Crank-Nicolson" };

  std::string str = "list of parameters:\n";
  str += "dim = " + d2s(DIM) + "\n";
  str += "scheme = " + time_scheme_name[TIME_SCHEME] + "\n";
  str += "mesh file name = " + MESH_FILE + "\n";
  //str += "mesh cl = " + d2s(CL) + "\n";
  str += "domain = [" + d2s(X_BEG) + ", " + d2s(X_END) + "] x [" + d2s(Y_BEG) + ", " + d2s(Y_END) + "]\n";
  str += "time = from " + d2s(TIME_BEG) + " to " + d2s(TIME_END) + " sec\n";
  str += "time step = " + d2s(TIME_STEP) + "\n";
  str += "number of time steps = " + d2s(N_TIME_STEPS) + "\n";
  str += "fe basis order = " + d2s(FE_ORDER) + "\n";
  str += "gamma = " + d2s(IPDG_GAMMA) + "\n";
  str += "print_vtu = " + d2s(PRINT_VTU) + "\n";
  str += "save_sol = " + d2s(SAVE_SOL) + "\n";
  str += "print_info = " + d2s(PRINT_INFO) + "\n";
  str += "vtu_step = " + d2s(VTU_STEP) + "\n";
  str += "sol_step = " + d2s(SOL_STEP) + "\n";
  str += "f0 = " + d2s(SOURCE_FREQUENCY) + "\n";
  str += "P = " + d2s(SOURCE_SUPPORT) + "\n";
  str += "xcen = " + d2s(SOURCE_CENTER_X) + "\n";
  str += "ycen = " + d2s(SOURCE_CENTER_Y) + "\n";

  str += "number of subdomains = " + d2s(N_SUBDOMAINS) + "\n";

  std::string coef_a = "", coef_b = "";
  for (int i = 0; i < N_SUBDOMAINS; ++i)
  {
    coef_a += (COEF_A_FILES[i] == "" ? d2s(COEF_A_VALUES[i]) : COEF_A_FILES[i]) + ", ";
    coef_b += (COEF_B_FILES[i] == "" ? d2s(COEF_B_VALUES[i]) : COEF_B_FILES[i]) + ", ";
  }

  str += "coef_a = " + coef_a + "\n";
  str += "coef_b = " + coef_b + "\n";
  return str;
}



void Parameters::establish_environment()
{
  generate_paths(); // generate all necessary path to files and directories
  check_clean_dirs(); // check the existance and clearance of some directories
}



void Parameters::generate_paths()
{
  MESH_FILE = MESH_DIR + "/" + MESH_FILE; // full path to the mesh file

  std::string coef_a = "", coef_b = "";
  for (int i = 0; i < N_SUBDOMAINS; ++i)
  {
    coef_a += (COEF_A_FILES[i] == "" ? d2s(COEF_A_VALUES[i]) : stem(COEF_A_FILES[i]));
    coef_b += (COEF_B_FILES[i] == "" ? d2s(COEF_B_VALUES[i]) : stem(COEF_B_FILES[i]));
  }

  RES_DIR = RES_TOP_DIR + "/" + stem(MESH_FILE) +
            "_T" + d2s(TIME_END) + "_K" + d2s(N_TIME_STEPS) +
            "_f" + d2s(SOURCE_FREQUENCY) + "_P" + d2s(SOURCE_SUPPORT) +
            "_A" + coef_a + "_B" + coef_b + "_gam" + d2s(IPDG_GAMMA) + "/";

  VTU_DIR = RES_DIR + "/" + VTU_DIR;
  SOL_DIR = RES_DIR + "/" + SOL_DIR;
  TIME_FILE = RES_DIR + "/" + TIME_FILE;
  INFO_FILE = RES_DIR + "/" + INFO_FILE;
}



void Parameters::check_clean_dirs() const
{
  using namespace boost::filesystem;

  path top_res_dir(RES_TOP_DIR); // top directory with results (from config.h, CMakeLists.txt)
  if (!exists(top_res_dir)) // if this top directory for results doesn't exist, we create it
    create_directory(top_res_dir);

  require(RES_DIR != "", "Directory for results has no name");
  path cur_res_dir(RES_DIR); // current directory with results
  if (exists(cur_res_dir) && is_directory(cur_res_dir)) // if this directory exists, we need to clean it up
    remove_all(cur_res_dir); // remove all contents of the directory and the directory itself
  create_directory(cur_res_dir); // now create empty directory

  if (PRINT_VTU)
    create_directory(VTU_DIR);

  if (SAVE_SOL)
    create_directory(SOL_DIR);
}



std::string stem(const std::string &filename)
{
  using namespace boost::filesystem;
  path filename_path(filename);
  return filename_path.stem().string();
}
