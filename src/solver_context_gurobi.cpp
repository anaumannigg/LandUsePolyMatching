#include "../include/solver_context.h"
#include "../include/solver_context_gurobi_access.h"
#include "gurobi_c++.h"

struct SolverContext::Impl {
    GRBEnv env;
    Impl() : env(true) {
        env.set("LogToConsole", "0");
        env.start();
    }
};

SolverContext::SolverContext() : impl_(std::make_unique<Impl>()) {}
SolverContext::~SolverContext() = default;

bool SolverContext::usesGurobi() const { return true; }

GRBEnv& SolverContextGurobiAccessor::env(SolverContext& ctx) {
    return ctx.impl_->env;
}
