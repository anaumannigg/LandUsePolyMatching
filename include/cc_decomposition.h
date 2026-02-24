#ifndef TCPOLYGONMATCHING_CC_DECOMPOSITION_H
#define TCPOLYGONMATCHING_CC_DECOMPOSITION_H
#include <vector>

#include "polygon_wh.h"

namespace CCDecomposition {
    void decomposeAndWriteToTXT(const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, const std::string& filename, const int& num_threads);

    void writeDecompositionToTXT(const int& num_computed_components, const std::vector<std::vector<int>>& comp1, const std::vector<std::vector<int>> comp2, const std::string& filename);
}

#endif //TCPOLYGONMATCHING_CC_DECOMPOSITION_H