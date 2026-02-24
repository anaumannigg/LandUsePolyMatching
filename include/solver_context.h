//
// Created by naumann on 22.01.26.
//

#ifndef TCPOLYGONMATCHING_SOLVER_CONTEXT_H
#define TCPOLYGONMATCHING_SOLVER_CONTEXT_H

#pragma once
#include <memory>

class SolverContext {
public:
    SolverContext();
    ~SolverContext();

    SolverContext(const SolverContext&) = delete;
    SolverContext& operator=(const SolverContext&) = delete;

    // optional: expose whether Gurobi is actually in use
    bool usesGurobi() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;

    friend struct SolverContextGurobiAccessor;
};

//singleton accessor
SolverContext& globalSolverContext();

#endif //TCPOLYGONMATCHING_SOLVER_CONTEXT_H