#include <angles/angles.hpp>

#include <vector>
#include <string>
#include <unordered_map>
#include <cctype>
#include <fstream>
#include <cstring>
#include <cmath>
#include <limits>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <fstream>

namespace angles {

Key toKey(std::string s) {
  for (char& ch : s ) ch = std::toupper(static_cast<unsigned char>(ch));
  if (s == "DES") return Key::DES;
  if (s == "DBS") return Key::DBS;
  if (s == "DEO") return Key::DEO;
  if (s == "DEA") return Key::DEA;
  return Key::Unknown;
}

Targets findTargets(const std::string& s, bool left) {
  Targets init{};
  switch ( toKey(s)) {
    case Key::DES:
      init = left ? Targets{1, 2, 19, 20, 21} : Targets{18, 17, 44, 43, 42}; break;
    case Key::DBS:
      init = left ? Targets{1, 2, 23, 24, 25} : Targets{22, 21, 56, 55, 54}; break;
    case Key::DEO:
      init = left ? Targets{1, 2, 11, 12, 13} : Targets{10, 9, 20, 19, 18}; break;
    case Key::DEA:
      init = left ? Targets{1, 2, 15, 16, 17} : Targets{14, 13, 32, 31, 30}; break;
    default:
      throw std::runtime_error("Unknown key");
  }

  init += -1;
  return init;
}


double massFromType(const std::vector<std::pair<int, double>>& type2mass, int type, double def ) {
  auto it = std::lower_bound(type2mass.begin(), type2mass.end(), type, [](const auto& p, int key){return p.first < key; });
  return (it != type2mass.end() && it->first == type) ? it->second : def;
}


double comZ(const std::vector<Atom>& atoms, const std::vector<int>& atomTypes, const std::vector<std::pair<int, double>>& type2mass){
  if (atoms.size() != atomTypes.size() ) {
    throw std::runtime_error("atoms/types size mismatch");
  }
  long double M = 0.0L, sz = 0.0L;
  for (size_t i=0; i<atoms.size(); i++) {
    const double m = massFromType(type2mass, atomTypes[i]);
    sz += m*atoms[i].z;
    M += m;
  }
  return (M == 0.0L) ? 0.0 : double(sz / M);
}

double wrapPBC(double x, double L) { // wrap x into [0, L)
  double y = std::fmod(x, L);
  if (y<0) y+=L;
  return y;
}

double comZPBC(const std::vector<Atom>& atoms, const std::vector<int>& atomTypes, const std::vector<std::pair<int, double>>& type2mass, double Lz){
  if (atoms.size() != atomTypes.size() ) {
    throw std::runtime_error("atoms/types size mismatch");
  }
  if(!(Lz > 0.0) ) {
    throw std::runtime_error("Lz must be positive");
  }
  long double C = 0.0L, S = 0.0L, M=0.0L;
  const long double twoPi = 2.0L * acosl(-1.0L);
  const long double k = twoPi / Lz;

  for (size_t i =0; i< atoms.size(); i++) {
    const double m = massFromType(type2mass, atomTypes[i]);
    const long double theta = k* wrapPBC(atoms[i].z, Lz);
    C += m * std::cos(theta);
    S += m * std::sin(theta);
    M += m;
  }
  if (M == 0.0L) return 0.0;

  long double ang = std::atan2(S, C);
  if (ang < 0) ang += twoPi;
  return double (ang * (Lz / twoPi));
}








// -------------- END STATISTICS ---------------

// -------------- START BASIC -------------------
double applyPBC(double x, double box){
  double hbox = box/2.0;
  double wrapped = fmod(x + hbox, box);
  if (wrapped < 0) wrapped += box;
  return wrapped - hbox;
}

double distance(const Atom& a, const Atom& b, const Box& box) {
  double dx, dy, dz, rsq;
  dx=applyPBC( a.x - b.x, box.x);
  dy=applyPBC( a.y - b.y, box.y);
  dz=applyPBC( a.z - b.z, box.z);
  rsq = dx*dx +dy*dy +dz*dz;
  return std::sqrt(rsq);
}

Angle computeTerminalAngle( Targets& target, const std::vector<Atom>& atoms, const Box& box) {
  Angle result;
  // compute C-C vector angle
  double dz;
  dz=applyPBC(atoms[target.t].z - atoms[target.c].z, box.z);
  result.ccCos = dz/distance(atoms[target.t], atoms[target.c], box);
  //std::cout << "target : " << target.t << " carbon : " << target.c << " dz " << dz << " result " << result.ccCos << "\n";

  // compute averaged C-H vector angle, assume CH bond length fixed
  double dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3, dx, dy;
  dx1 = applyPBC(atoms[target.h1].x - atoms[target.t].x, box.x);
  dx2 = applyPBC(atoms[target.h2].x - atoms[target.t].x, box.x);
  dx3 = applyPBC(atoms[target.h3].x - atoms[target.t].x, box.x);
  dx = dx1+dx2+dx3;

  dy1 = applyPBC(atoms[target.h1].y - atoms[target.t].y, box.y);
  dy2 = applyPBC(atoms[target.h2].y - atoms[target.t].y, box.y);
  dy3 = applyPBC(atoms[target.h3].y - atoms[target.t].y, box.y);
  dy = dy1 + dy2 + dy3;

  dz1 = applyPBC(atoms[target.h1].z - atoms[target.t].z, box.z);
  dz2 = applyPBC(atoms[target.h2].z - atoms[target.t].z, box.z);
  dz3 = applyPBC(atoms[target.h3].z - atoms[target.t].z, box.z);
  dz = dz1 + dz2 + dz3;

  double rsq = dx*dx + dy*dy + dz*dz;
  result.chCos = dz/std::sqrt(rsq);

  result.zloc = atoms[target.t].z;
  
  return result;
}

Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y*b.z - a.z*b.y,
          a.z*b.x - a.x*b.z,
          a.x*b.y - a.y*b.x};
}
double dot(const Vec3& a, const Vec3& b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}
double norm(const Vec3& a) {
  return std::sqrt(dot(a,a));
}

Dihedral computeDihedralAngle( DTargets& target, const std::vector<Atom>& atoms, const Box& box) {
  const Vec3& r1{atoms[target.d1].x, atoms[target.d1].y, atoms[target.d1].z};
  const Vec3& r2{atoms[target.d2].x, atoms[target.d2].y, atoms[target.d2].z};
  const Vec3& r3{atoms[target.d3].x, atoms[target.d3].y, atoms[target.d3].z};
  const Vec3& r4{atoms[target.d4].x, atoms[target.d4].y, atoms[target.d4].z};

  Dihedral out;
  // mean z-location
  out.zloc = 0.25* ( r1.z + r2.z + r3.z + r4.z);

  // bond vectors
  Vec3 b1{applyPBC(r2.x-r1.x, box.x), applyPBC(r2.y-r1.y, box.y), applyPBC(r2.z-r1.z, box.z)};
  Vec3 b2{applyPBC(r3.x-r2.x, box.x), applyPBC(r3.y-r2.y, box.y), applyPBC(r3.z-r2.z, box.z)};
  Vec3 b3{applyPBC(r4.x-r3.x, box.x), applyPBC(r4.y-r3.y, box.y), applyPBC(r4.z-r3.z, box.z)};

  // plane normal
  Vec3 n1=cross(b1, b2);
  Vec3 n2=cross(b2, b3);

  // angle from atan2
  double b2n = norm(b2);
  if (b2n==0.0) {
    out.angle=0.0;
    return out;
  }
  Vec3 b2hat{ b2.x/b2n, b2.y/b2n, b2.z/b2n};
  double x = dot(n1, n2);
  double y = dot(cross(n1, n2), b2hat);
  out.angle = std::atan2(y, x); // -pi to pi
  return out;
}

} // -----------namespace angles




