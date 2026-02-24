#include "../include/solver_context.h"

SolverContext& globalSolverContext() {
    static SolverContext ctx;
    return ctx;
}
