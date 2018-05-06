#include "LHEreader.h"

namespace myFuncs {
LHEreader::LHEreader(const std::string& filename) {
  m_instream = std::ifstream(filename);
  if (!m_instream.is_open()) {
    throw "Can't open file " + filename;
  }
}

std::vector<Particle> LHEreader::getEventParticles() {
  std::vector<Particle> particles;
  particles.reserve(6);

  std::string line;
  // Seek <event>
  while (!m_instream.eof()) {

    std::getline(m_instream, line);

    // Continue until <event>
    if (line.find("<event>") != std::string::npos) {
      break;
    }
  }

  // ignore event info line
  std::getline(m_instream, line);

  // Found <event>
  while (true) {
    Particle particle;
    int iDummy;
    double dDummy;
    m_instream >> particle.pdg >> iDummy >> iDummy >> iDummy >> iDummy >> iDummy >> particle.px >> particle.py >> particle.pz >>
        particle.E >> particle.mass >> dDummy >> dDummy;

    // Couldn't parse - either </event> or bad line
    if (m_instream.fail()) {
      m_instream.clear();
      return particles;
    }

    // Add particle
    particles.push_back(particle);

    // Seek line end
    m_instream.ignore(1000, '\n');
  }

  return particles;
}

} // namespace myFuncs
