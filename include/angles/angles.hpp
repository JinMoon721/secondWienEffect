#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>

#include <ioTraj/atom.hpp>
using namespace ioTraj;

// ----------------- namespace angles -------------------------------
namespace angles {


struct Angle { 
  double zloc;  // z location of terminal carbon
  double ccCos; // angle from +z axis of terminal cc bond
  double chCos; // angle from +z axis of the mean terminal ch bonds
};


struct Dihedral { 
  double zloc;  // z location of dihedral pairs
  double angle; // radian angle from -pi to pi
};

struct Targets {
  int t; // terminal 
  int c; // carbon next to terminal
  int h1; // hydrogen1 at terminal
  int h2;
  int h3;
  Targets operator+(int k) const { return {t+k, c+k, h1+k, h2+k, h3+k}; }
  Targets& operator+=(int k) { t+=k; c+=k; h1+=k; h2+=k; h3+=k; return *this; }
};

struct DTargets { // dihedral target
  int d1; 
  int d2; 
  int d3; 
  int d4;
  DTargets operator+(int k) const { return {d1+k, d2+k, d3+k, d4+k}; }
  DTargets& operator+=(int k) { d1+=k; d2+=k; d3+=k; d4+=k; return *this; }
};

enum class Key { DES, DBS, DEO, DEA, Unknown };
Key toKey(std::string s); 

Targets findTargets(const std::string& s, bool left); 
double massFromType(const std::vector<std::pair<int, double>>& type2mass, int type, double def = -1.0);
double comZ(const std::vector<Atom>& atoms, const std::vector<int>& atomTypes, const std::vector<std::pair<int, double>>& type2mass);
double comZPBC(const std::vector<Atom>& atoms, const std::vector<int>& atomTypes, const std::vector<std::pair<int, double>>& type2mass, double Lz);

Angle computeTerminalAngle( Targets& target, const std::vector<Atom>& atoms, const Box& box);

struct Vec3 {
  double x, y, z;
};
Vec3 cross(const Vec3& a, const Vec3& b);
double dot(const Vec3& a, const Vec3& b);
double norm(const Vec3& a);
Dihedral computeDihedralAngle( DTargets& target, const std::vector<Atom>& atoms, const Box& box);

double applyPBC(double x, double box);
double distance(const Atom& a, const Atom& b, const Box& box);
double wrapPBC(double x, double L);

} // -----------namespace angles


