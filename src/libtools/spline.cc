#include <votca/tools/spline.h>

namespace votca {
    namespace tools {

        using namespace std;

        int Spline::GenerateGrid(double min, double max, double h) {
            int vec_size = (int) ((max - min) / h + 1.00000001);
            _r.resize(vec_size);
            int i;

            double r_init;

            for (r_init = min, i = 0; i < vec_size - 1; r_init += h) {
                _r[i++] = r_init;
            }
            _r[i] = max;
            _f.resize(_r.size(), false);
            _f2.resize(_r.size(), false);
            return _r.size();
        }

    }
}
