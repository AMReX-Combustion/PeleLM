
#include <pelelm_prob.H>

namespace ProbParm
{
   AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean = 101325.0;
   AMREX_GPU_DEVICE_MANAGED  amrex::Real standoff = - 0.022;
   AMREX_GPU_DEVICE_MANAGED  amrex::Real pertmag = 0.0004;

   AMREX_GPU_DEVICE_MANAGED unsigned int pmf_N = 0;
   AMREX_GPU_DEVICE_MANAGED unsigned int pmf_M = 0;
   AMREX_GPU_DEVICE_MANAGED bool pmf_do_average = false;

   amrex::Gpu::ManagedVector<amrex::Real>* pmf_X = nullptr;
   amrex::Gpu::ManagedVector<amrex::Real>* pmf_Y = nullptr;
   amrex::Gpu::ManagedVector<amrex::Real>* fuel_state = nullptr;

   AMREX_GPU_DEVICE_MANAGED amrex::Real* d_pmf_X = nullptr;
   AMREX_GPU_DEVICE_MANAGED amrex::Real* d_pmf_Y = nullptr;
   AMREX_GPU_DEVICE_MANAGED amrex::Real* d_fuel_state = nullptr;

   std::string pmf_datafile = "";
   amrex::Vector<std::string> pmf_names;
} // namespace ProbParm

std::string
read_pmf_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

bool
checkQuotes(const std::string& str)
{
  int count = 0;
  for (char c : str) {
    if (c == '"')
      count++;
  }
  if ((count % 2) == 0)
    return true;
  else
    return false;
}

void
read_pmf(const std::string myfile)
{
  std::string firstline, secondline, remaininglines;
  int pos1, pos2;
  int variable_count, line_count;

  std::ifstream infile(myfile);
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  if (!checkQuotes(firstline))
    amrex::Abort("PMF file variable quotes unbalanced");
  std::getline(iss, secondline);
  pos1 = 0;
  pos2 = 0;
  variable_count = 0;
  while ((pos1 < firstline.length() - 1) && (pos2 < firstline.length() - 1)) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    variable_count++;
    pos1 = pos2 + 1;
  }

  ProbParm::pmf_names.resize(variable_count);
  pos1 = 0;
  pos2 = 0;
  for (int i = 0; i < variable_count; i++) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    ProbParm::pmf_names[i] = firstline.substr(pos1 + 1, pos2 - (pos1 + 1));
    pos1 = pos2 + 1;
  }

  amrex::Print() << variable_count << " variables found in PMF file"
                 << std::endl;
  // for (int i = 0; i < variable_count; i++)
  //  amrex::Print() << "Variable found: " << ProbParm::pmf_names[i] <<
  //  std::endl;

  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  ProbParm::pmf_N = line_count;
  ProbParm::pmf_M = variable_count - 1;
  ProbParm::pmf_X->resize(ProbParm::pmf_N);
  ProbParm::pmf_Y->resize(ProbParm::pmf_N * ProbParm::pmf_M);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (int i = 0; i < ProbParm::pmf_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> (*ProbParm::pmf_X)[i];
    for (int j = 0; j < ProbParm::pmf_M; j++) {
      sinput >> (*ProbParm::pmf_Y)[j * ProbParm::pmf_N + i];
    }
  }
  ProbParm::d_pmf_X = ProbParm::pmf_X->dataPtr();
  ProbParm::d_pmf_Y = ProbParm::pmf_Y->dataPtr();
}

extern "C" {
    void amrex_probinit (const int* init,
                         const int* name,
                         const int* namelen,
                         const amrex_real* problo,
                         const amrex_real* probhi)
    {
        amrex::ParmParse pp("prob");

        pp.query("P_mean", ProbParm::P_mean);
        pp.query("standoff", ProbParm::standoff);
        pp.query("pertmag", ProbParm::pertmag);

        pp.query("pmf_datafile", ProbParm::pmf_datafile);

        ProbParm::pmf_X = new amrex::Gpu::ManagedVector<amrex::Real>;
        ProbParm::pmf_Y = new amrex::Gpu::ManagedVector<amrex::Real>;

        read_pmf(ProbParm::pmf_datafile);

    }
}
