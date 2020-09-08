#ifndef PRESSURE_H_
#define PRESSURE_H_

#include "result.h"

namespace Results {
    class Pressure: public Result
    {
    public:
        Pressure(const std::string& filename);
        void notify(int current_step, double current_distance, const SimulationData& data) override;
        void finalize() override;

    private:
        std::string filename;
    };
}

#endif
