#include <iostream>
#include <stdio.h>
#include <dolfin.h>
#include "stokes.h"
#include "fasp_solver.h"
using namespace FASPFSI;
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
	int M = 32;
	auto mesh = std::make_shared<UnitSquareMesh>(M, M);
	auto W = std::make_shared<stokes::FunctionSpace>(mesh);
	parameters["linear_algebra_backend"] = "Eigen";	

	double dt = 0.002;
	double T = 0.1;

	auto u0 = std::make_shared<Constant>(0.0, 0.0);
	auto u1 = std::make_shared<Constant>(1.0, 0.0);
	auto boundary0 = std::make_shared<LRB_Boundary>();
	auto boundary1 = std::make_shared<T_Boundary>();
//	auto bc0 = make_shared<DirichletBC>(W->sub(0), u0, boundary0);	
//	auto bc1 = make_shared<DirichletBC>(W->sub(0), u1, boundary1);
	DirichletBC bc0(W->sub(0), u0, boundary0);	
	DirichletBC bc1(W->sub(0), u1, boundary1);
	std::vector<const DirichletBC*> bcs{{&bc0, &bc1}};	
	auto f = std::make_shared<Constant>(0.0, 0.0);
	stokes::BilinearForm a(W, W);
	stokes::LinearForm L(W);
	L.f = f;
	auto inv_dt = std::make_shared<Constant>(1.0/dt);
	a.inv_dt = inv_dt;
	L.inv_dt = inv_dt;
	auto w = std::make_shared<Function>(W);	
	auto w0 = std::make_shared<Function>(W);
	L.w0 = w0;
	FaspSolver faspsolver;
	auto A = std::make_shared<EigenMatrix>();
	EigenVector b;
	EigenVector *x = (EigenVector*)&(*w->vector());
	EigenVector *x0 = (EigenVector*)&(*w0->vector());
	double *xval = (double *)x->data();
	double *x0val = (double *)x0->data();
	assemble(*A, a);
	
	faspsolver.assign_idx(&(*W));
	double t = dt;
	File ufile_pvd("result/velocity.pvd");
	File pfile_pvd("result/pressure.pvd");

	while (t < T + DOLFIN_EPS)
	{
		assemble(b, L);
		bcs[0]->apply(*A, b);
		bcs[1]->apply(*A, b);
		if(t == dt){
			faspsolver.set_operator(*A);
		}
		faspsolver.solve(*x, b, 1, 0);
	
		ufile_pvd << (*w)[0];
		pfile_pvd << (*w)[1];
		for(int i = 0; i < x->size(); ++i)
		{
			x0val[i] = xval[i];
		}
		t+=dt;		
	}
	plot((*w)[0]);
	interactive();
	return 0;	
}
