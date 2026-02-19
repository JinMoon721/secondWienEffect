#include <ioTraj/atom.hpp>
#include <angles/angles.hpp>
#include <ioTraj/dcdreader.hpp>
#include <ioTraj/lammpsreader.hpp>
#include <stats/stats.hpp>

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace angles;
using namespace ioTraj;
using namespace stats;

void user_printBinary(std::string& file,  std::vector<std::vector<float>>& matrix){
  std::cout << "Binary Output File generated: " << file << std::endl;
  std::cout << "Rows: " << matrix.size() << ", Columns: " << (matrix.empty() ? 0 : matrix[0].size()) << std::endl;

  std::ofstream out ( file, std::ios::binary);
  int rows = matrix.size();
  int cols = matrix[0].size();
  out.write(reinterpret_cast<const char*>(&rows), sizeof(int));
  out.write(reinterpret_cast<const char*>(&cols), sizeof(int));

  for (const auto& row : matrix) {
    out.write(reinterpret_cast<const char*>(row.data()), sizeof(float) * cols);
  }

  out.close();
}

int main(int argc, char** argv) {
  if (argc != 4 ) {
    std::cerr << "Error: Num Args not Match\n" ;
    std::cerr << "Usage : " << argv[0] << " dirName density field \n";
    return 1;
  }
  std::string dirName=argv[1];
  std::string density= argv[2];
  std::string field = argv[3];

  // Read DCD file
  std::string dump = std::string("../data/dumps/") + dirName + "/dumpD"  + density + "E" + field + ".dcd";
  std::cout << "---- Reading dcd file at " << dump << " ----\n";
  DCDReader reader(dump);
  int numatoms = reader.natoms();
  std::cout << "number of atoms: " << numatoms << "\n";
  std::cout << "number of frames (original): " << reader.nframes() << "\n";
  auto frames = reader.read_all();
  size_t numsnap = frames.size();
  std::cout << "Read Traj of " << numsnap <<  " frames\n";

  std::vector<std::vector<float>> out;
  std::string outf=std::string("../data/traj/") + dirName + "/trajD" +density + "E" + field + ".binary"; 

  for (size_t time=0; time<numsnap; time++ ) { 
    auto frame = frames[time].atoms;
    for (int atom=0; atom < numatoms; atom++ ) {
      float type; // original dcd file have ordered 
      if ( atom < numatoms/2) {
        type = 3;
      } else {
        type =4;
      }
      std::vector<float> one;
      one.push_back(type);
      one.push_back( frame[atom].x);
      one.push_back( frame[atom].y);
      one.push_back( frame[atom].z);
      out.push_back(std::move(one));
    }
  }
  user_printBinary(outf, out);
  return 0;
}

