#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace myFuncs {

struct Particle {
  int pdg;
  double px;
  double py;
  double pz;
  double E;
  double mass;
};

std::ostream& operator<<(std::ostream& ostream, Particle const& particle) {
  return ostream << particle.pdg << "\t" << particle.px << "\t" << particle.py << "\t" << particle.pz << "\t" << particle.E << "\t"
          << particle.mass << "\n";
}

class LHEreader {
public:
  // Constructor
  LHEreader(const std::string& filename);

  // Read particles in next event
  std::vector<Particle> getEventParticles();

private:
  std::ifstream m_instream;
};

} // namespace myFuncs
