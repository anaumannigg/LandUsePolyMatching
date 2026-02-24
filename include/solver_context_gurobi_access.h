#ifndef TCPOLYGONMATCHING_SOLVER_CONTEXT_GUROBI_ACCESS_H
#define TCPOLYGONMATCHING_SOLVER_CONTEXT_GUROBI_ACCESS_H

#pragma once
class SolverContext;
class GRBEnv;

struct SolverContextGurobiAccessor {
    static GRBEnv& env(SolverContext& ctx);
};


#endif //TCPOLYGONMATCHING_SOLVER_CONTEXT_GUROBI_ACCESS_H