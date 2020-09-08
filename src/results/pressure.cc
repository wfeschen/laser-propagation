#include "pressure.h"
#include "../util/util.h"
#include "../util/io.h"

namespace Results {
  Pressure::Pressure(const std::string& filename)
    :filename(filename) {
    IO::clear_contents(filename);
  }
 
  void Pressure::notify(int, double distance, const SimulationData& data) {
    IO::write_append(filename, distance, data.pressure);
  }
 
  void Pressure::finalize() {}
}
