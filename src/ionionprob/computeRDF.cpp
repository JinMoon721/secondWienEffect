#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <iomanip>
#include <map>

#include <stats/stats.hpp>
using namespace stats;




std::vector<std::vector<double>> readBTrajFromFile(const std::string& filename) {
  std::ifstream in(filename, std::ios::binary);
  int32_t rows32=0, cols32=0;
  in.read(reinterpret_cast<char*>(&rows32), sizeof(rows32));
  in.read(reinterpret_cast<char*>(&cols32), sizeof(cols32));

  size_t rows = static_cast<size_t>(rows32);
  size_t cols = static_cast<size_t>(cols32);

  std::vector<std::vector<double>> out(rows, std::vector<double>(cols));
  std::vector<float> rowf(cols);
  for( size_t i=0; i<rows; i++){
    in.read(reinterpret_cast<char*>(rowf.data()), static_cast<std::streamsize>( sizeof(float)*cols));
    for (size_t j=0; j<cols; j++) {
      out[i][j] = static_cast<double>(rowf[j]);
    }
   }
  in.close();
  std::cout << "file reading done" << std::endl;
  return out;
}







float fieldConvert(const std::string& s) {
  if (s.find_first_not_of('0') == std::string::npos) {
    return 0.0f;
  }
  size_t leadingZeros = s.find_first_not_of('0');

  if ( leadingZeros >= 2) {
    std::string digits = s.substr(leadingZeros);
    float value = std::stof(digits);
    return value / std::pow(10, leadingZeros -1);
  }
  return std::stof(s);
}

int main(int argc, char* argv[]) {
  if (argc != 4 ) {
    std::cerr << "Error: Not enough arguments.\n" ;
    std::cerr << "Usage : " << argv[0] << " dir_name density field\n";
    return 1;
  }
  std::string dirName=argv[1]; // name of directory to output results, in results directory
  std::string rho = argv[2]; 
  int n = rho.size();
  long long val = std::stoll(rho);
  float density = static_cast<float>(val) / std::pow(10.0, n-1);
  std::string fieldname = argv[3];
  double field = fieldConvert(fieldname)* 0.023; //kcal /molA,
  double fconversion = 25.7/0.592 ; // kcal/molA to mV/A

  // get dist data, flatten it to feed it to paralle job
  std::string distFile = std::string("../data/cnnDist/") + dirName + "D" + rho + "E" + fieldname + ".binary";
  auto dist=readBTrajFromFile(distFile);
  size_t numsnap = dist.size();
  size_t numatoms = dist[0].size();

  std::cout << "\nAnalysis Setups\n";
  std::cout << "data location: " << distFile << std::endl ; 
  std::cout << "numIons: " << numatoms << " numPairs: " << numatoms/2 << "\n";
  std::cout << "density of ions : " << density << "M \n";
  std::cout << "fieldStrength: " << field << " kcal/molA = " << field * fconversion << " mV/A" << "\n"; 
  std::cout << "\n\n";

  std::cout << "Radial Distribution Function computation Starts...\n";

  //flatten 2d double std vector into 1d vector
  std::vector<double> fdist;
  for (auto& row : dist) {
    fdist.insert(fdist.end(), row.begin(), row.end());
  }

  /// compute radial histograms
  size_t numBins=1000;
  auto rdf = makeHistogramAuto(fdist, numBins);

  std::string outfile = std::string("../results/ionionprob/") + dirName + "D" + rho + "E" +  fieldname + ".dat";
  std::ofstream out(outfile);
  writeHistogram(rdf, out);

  std::cout << "Radial Distribution Function computation Ended.\n";
  return 0;
}

  



