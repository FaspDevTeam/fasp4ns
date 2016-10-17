#include <iostream>
#include <stdio.h>
#include <dolfin.h>
#include "stokes.h"
#include "fasp_solver.h"
using namespace FASPFSI;
using namespace std;
using namespace dolfin;
class LRB_Boundary : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS
			or x[1] < DOLFIN_EPS;
	}
};
class T_Boundary : public SubDomain
{
	bool inside(const Array<double>& x, bool on_boundary) const
	{
		return x[1] > 1 - DOLFIN_EPS;
	}
};

int main(int argc, char* argv[]){
	auto mesh = make_shared<UnitSquareMesh>(32, 32);
	auto W = make_shared<stokes::FunctionSpace>(mesh);
	parameters["linear_algebra_backend"] = "Eigen";	
	auto u0 = make_shared<Constant>(0.0, 0.0);
	auto u1 = make_shared<Constant>(1.0, 0.0);
	auto boundary0 = make_shared<LRB_Boundary>();
	auto boundary1 = make_shared<T_Boundary>();
//	auto bc0 = make_shared<DirichletBC>(W->sub(0), u0, boundary0);	
//	auto bc1 = make_shared<DirichletBC>(W->sub(0), u1, boundary1);
	DirichletBC bc0(W->sub(0), u0, boundary0);	
	DirichletBC bc1(W->sub(0), u1, boundary1);
	vector<const DirichletBC*> bcs{{&bc0, &bc1}};	
	auto f = make_shared<Constant>(0.0, 0.0);
	stokes::BilinearForm a(W, W);
	stokes::LinearForm L(W);
	L.f = f;
	
	auto w = std::make_shared<Function>(W);	
	FaspSolver faspsolver;
	auto A = make_shared<EigenMatrix>();
	EigenVector b;
	EigenVector *x = (EigenVector*)&(*w->vector());
	assemble(*A, a);
	assemble(b, L);
	bcs[0]->apply(*A, b);
	bcs[1]->apply(*A, b);
//	bcs[0]->apply(*w->vector());
//	bcs[1]->apply(*w->vector());
//	solve(a == L, w, bcs);
//	solve(*A, *w->vector(), b);
	faspsolver.assign_idx(&(*W));
//	faspsolver.assemble_blk();
	faspsolver.set_operator(*A);
	faspsolver.solve(*x, b, 1, 1);
//	faspsolver.print_FS();
	File ufile_pvd("velocity.pvd");
	File pfile_pvd("pressure.pvd");
	Function u = (*w)[0];
	Function p = (*w)[1];
	ufile_pvd << u;
	pfile_pvd << p;

	plot(u);
	interactive();
	return 0;	
}
