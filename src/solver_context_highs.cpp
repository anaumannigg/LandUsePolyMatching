#include "../include/solver_context.h"

struct SolverContext::Impl {
    // empty, context is only needed for Gurobi, this is a placeholder
};

SolverContext::SolverContext() : impl_(std::make_unique<Impl>()) {}
SolverContext::~SolverContext() = default;

bool SolverContext::usesGurobi() const { return false; }
