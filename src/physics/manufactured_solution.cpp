#include <Sacado.hpp>
#include <deal.II/base/function.h>
#include <deal.II/base/function.templates.h> // Needed to instantiate dealii::Function<PHILIP_DIM,Sacado::Fad::DFad<double>>
#include <deal.II/base/function_time.templates.h> // Needed to instantiate dealii::Function<PHILIP_DIM,Sacado::Fad::DFad<double>>

#include "manufactured_solution.h"

template class dealii::FunctionTime<Sacado::Fad::DFad<double>>; // Needed by Function
template class dealii::Function<PHILIP_DIM,Sacado::Fad::DFad<double>>;

namespace PHiLiP {

bool isfinite(Sacado::Fad::DFad<double> value)
{
    return std::isfinite(static_cast<double>(value.val()));
}

template <int dim, typename real>
inline real ManufacturedSolutionSine<dim,real>
::value (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    real value = this->amplitudes[istate];
    for (int d=0; d<dim; d++) {
        value *= sin( this->frequencies[istate][d] * point[d] );
        assert(isfinite(value));
    }
    value += this->base_values[istate];
    return value;
}

template <int dim, typename real>
inline real ManufacturedSolutionAdd<dim,real>
::value (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    real value = 0.0;
    for (int d=0; d<dim; d++) {
        value += this->amplitudes[istate]*sin( this->frequencies[istate][d] * point[d] );
        assert(isfinite(value));
    }
    value += this->base_values[istate];
    return value;
}

template <int dim, typename real>
inline real ManufacturedSolutionCosine<dim,real>
::value (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    real value = this->amplitudes[istate];
    for (int d=0; d<dim; d++) {
        value *= cos( this->frequencies[istate][d] * point[d] );
        assert(isfinite(value));
    }
    value += this->base_values[istate];
    return value;
}

template <int dim, typename real>
inline real ManufacturedSolutionExp<dim,real>
::value (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    real value = 0.0;
    for (int d=0; d<dim; d++) {
        value += exp( point[d] );
        assert(isfinite(value));
    }
    value += this->base_values[istate];
    return value;
}

template <int dim, typename real>
inline real ManufacturedSolutionEvenPoly<dim,real>
::value (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    real value = 0.0;
    const double poly_max = 7;
    for (int d=0; d<dim; d++) {
        value += pow(point[d] + 0.5, poly_max);
    }
    value += this->base_values[istate];
    return value;
}

template <int dim, typename real>
inline real ManufacturedSolutionPoly<dim,real>
::value (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    real value = 0.0;
    for (int d=0; d<dim; d++) {
        const real x = point[d];
        value += 1.0 + x - x*x - x*x*x + x*x*x*x - x*x*x*x*x + x*x*x*x*x*x + 0.001*sin(50*x);
    }
    value += this->base_values[istate];
    return value;
}

template <int dim, typename real>
inline real ManufacturedSolutionAtan<dim,real>
::value(const dealii::Point<dim,real> &point, const unsigned int /*istate*/) const
{
    real val = 1.0;
    for(unsigned int i = 0; i < dim; ++i){
        real x = point[i];
        real val_dim = 0;
        for(unsigned int j = 0; j < n_shocks[i]; ++j){
            // taking the product of function in each direction
            val_dim += atan(S_j[i][j]*(x-x_j[i][j]));
        }
        val *= val_dim;
    }
    return val;
}

template <int dim, typename real>
inline real ManufacturedSolutionBoundaryLayer<dim,real>
::value(const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    real val = 1.0;
    for(unsigned int d = 0; d < dim; ++d){
        real x = point[d];
        val *= x + (exp(x/epsilon[istate][d])-1.0)/(1.0-exp(1.0/epsilon[istate][d]));
    }
    return val;
}

template <int dim, typename real>
inline dealii::Tensor<1,dim,real> ManufacturedSolutionSine<dim,real>
::gradient (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    dealii::Tensor<1,dim,real> gradient;
    for (int dim_deri=0; dim_deri<dim; dim_deri++) {
        gradient[dim_deri] = this->amplitudes[istate] * this->frequencies[istate][dim_deri];
        for (int dim_trig=0; dim_trig<dim; dim_trig++) {
            const real angle = this->frequencies[istate][dim_trig] * point[dim_trig];
            if (dim_deri == dim_trig) gradient[dim_deri] *= cos( angle );
            if (dim_deri != dim_trig) gradient[dim_deri] *= sin( angle );
        }
        assert(isfinite(gradient[dim_deri]));
    }
    // Hard-coded is much more readable than the dimensionally generic one
    const real A = this->amplitudes[istate];
    const dealii::Tensor<1,dim,real> f = this->frequencies[istate];
    if (dim==1) {
        const real fx = f[0]*point[0];
        gradient[0] = A*f[0]*cos(fx);
    }
    if (dim==2) {
        const real fx = f[0]*point[0];
        const real fy = f[1]*point[1];
        gradient[0] = A*f[0]*cos(fx)*sin(fy);
        gradient[1] = A*f[1]*sin(fx)*cos(fy);
    }
    if (dim==3) {
        const real fx = f[0]*point[0];
        const real fy = f[1]*point[1];
        const real fz = f[2]*point[2];
        gradient[0] = A*f[0]*cos(fx)*sin(fy)*sin(fz);
        gradient[1] = A*f[1]*sin(fx)*cos(fy)*sin(fz);
        gradient[2] = A*f[2]*sin(fx)*sin(fy)*cos(fz);
    }
    return gradient;
}

template <int dim, typename real>
inline dealii::Tensor<1,dim,real> ManufacturedSolutionAdd<dim,real>
::gradient (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    dealii::Tensor<1,dim,real> gradient;
    const real A = this->amplitudes[istate];
    const dealii::Tensor<1,dim,real> f = this->frequencies[istate];
    if (dim==1) {
        const real fx = f[0]*point[0];
        gradient[0] = A*f[0]*cos(fx);
    }
    if (dim==2) {
        const real fx = f[0]*point[0];
        const real fy = f[1]*point[1];
        gradient[0] = A*f[0]*cos(fx);
        gradient[1] = A*f[1]*cos(fy);
    }
    if (dim==3) {
        const real fx = f[0]*point[0];
        const real fy = f[1]*point[1];
        const real fz = f[2]*point[2];
        gradient[0] = A*f[0]*cos(fx);
        gradient[1] = A*f[1]*cos(fy);
        gradient[2] = A*f[2]*cos(fz);
    }
    return gradient;
}

template <int dim, typename real>
inline dealii::Tensor<1,dim,real> ManufacturedSolutionCosine<dim,real>
::gradient (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    dealii::Tensor<1,dim,real> gradient;
    const real A = this->amplitudes[istate];
    const dealii::Tensor<1,dim,real> f = this->frequencies[istate];
    if (dim==1) {
        const real fx = f[0]*point[0];
        gradient[0] = -A*f[0]*sin(fx);
    }
    if (dim==2) {
        const real fx = f[0]*point[0];
        const real fy = f[1]*point[1];
        gradient[0] = -A*f[0]*sin(fx)*cos(fy);
        gradient[1] = -A*f[1]*cos(fx)*sin(fy);
    }
    if (dim==3) {
        const real fx = f[0]*point[0];
        const real fy = f[1]*point[1];
        const real fz = f[2]*point[2];
        gradient[0] = -A*f[0]*sin(fx)*cos(fy)*cos(fz);
        gradient[1] = -A*f[1]*cos(fx)*sin(fy)*cos(fz);
        gradient[2] = -A*f[2]*cos(fx)*cos(fy)*sin(fz);
    }
    return gradient;
}

template <int dim, typename real>
inline dealii::Tensor<1,dim,real> ManufacturedSolutionExp<dim,real>
::gradient (const dealii::Point<dim,real> &point, const unsigned int /*istate*/) const
{
    dealii::Tensor<1,dim,real> gradient;
    if (dim==1) {
        gradient[0] = exp(point[0]);
    }
    if (dim==2) {
        gradient[0] = exp(point[0]);
        gradient[1] = exp(point[1]);
    }
    if (dim==3) {
        gradient[0] = exp(point[0]);
        gradient[1] = exp(point[1]);
        gradient[2] = exp(point[2]);
    }
    return gradient;
}

template <int dim, typename real>
inline dealii::Tensor<1,dim,real> ManufacturedSolutionEvenPoly<dim,real>
::gradient (const dealii::Point<dim,real> &point, const unsigned int  /*istate*/) const
{
    dealii::Tensor<1,dim,real> gradient;
    const double poly_max = 7;
    if (dim==1) {
        gradient[0] = poly_max*pow(point[0] + 0.5, poly_max-1);
    }
    if (dim==2) {
        gradient[0] = poly_max*pow(point[0] + 0.5, poly_max-1);
        gradient[1] = poly_max*pow(point[1] + 0.5, poly_max-1);
    }
    if (dim==3) {
        gradient[0] = poly_max*pow(point[0] + 0.5, poly_max-1);
        gradient[1] = poly_max*pow(point[1] + 0.5, poly_max-1);
        gradient[2] = poly_max*pow(point[2] + 0.5, poly_max-1);
    }
    return gradient;
}

template <int dim, typename real>
inline dealii::Tensor<1,dim,real> ManufacturedSolutionPoly<dim,real>
::gradient (const dealii::Point<dim,real> &point, const unsigned int  /*istate*/) const
{
    dealii::Tensor<1,dim,real> gradient;
    if (dim==1) {
        const real x = point[0];
        gradient[0] = 1.0 - 2*x -3*x*x + 4*x*x*x - 5*x*x*x*x + 6*x*x*x*x*x + 0.050*cos(50*x);
    }
    if (dim==2) {
        real x = point[0];
        gradient[0] = 1.0 - 2*x -3*x*x + 4*x*x*x - 5*x*x*x*x + 6*x*x*x*x*x + 0.050*cos(50*x);
        x = point[1];
        gradient[1] = 1.0 - 2*x -3*x*x + 4*x*x*x - 5*x*x*x*x + 6*x*x*x*x*x + 0.050*cos(50*x);
    }
    if (dim==3) {
        real x = point[0];
        gradient[0] = 1.0 - 2*x -3*x*x + 4*x*x*x - 5*x*x*x*x + 6*x*x*x*x*x;
        x = point[1];
        gradient[1] = 1.0 - 2*x -3*x*x + 4*x*x*x - 5*x*x*x*x + 6*x*x*x*x*x;
        x = point[2];
        gradient[2] = 1.0 - 2*x -3*x*x + 4*x*x*x - 5*x*x*x*x + 6*x*x*x*x*x;
    }
    return gradient;
}

template <int dim, typename real>
inline dealii::Tensor<1,dim,real> ManufacturedSolutionAtan<dim,real>
::gradient(const dealii::Point<dim,real> &point, const unsigned int /*istate*/) const
{
    dealii::Tensor<1,dim,real> gradient;
    for(unsigned int k = 0; k < dim; ++k){
        // taking the k^th derivative
        real grad_dim = 1;
        for(unsigned int i = 0; i < dim; ++i){
            real x = point[i];
            real val_dim = 0;
            for(unsigned int j = 0; j < n_shocks[i]; ++j){
                if(i==k){
                    // taking the derivative dimension
                    real coeff = S_j[i][j]*(x-x_j[i][j]);
                    val_dim += S_j[i][j]/(pow(coeff,2)+1);
                }else{
                    // value product unaffected
                    val_dim += atan(S_j[i][j]*(x-x_j[i][j]));
                }
            }
            grad_dim *= val_dim;
        }
        gradient[k] = grad_dim;
    }
    return gradient;
}

template <int dim, typename real>
inline dealii::Tensor<1,dim,real> ManufacturedSolutionBoundaryLayer<dim,real>
::gradient(const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    dealii::Tensor<1,dim,real> gradient;
    if(dim == 1){
        const real x = point[0];
        gradient[0] = (1 + (exp(x/epsilon[istate][0])/epsilon[istate][0])/(1.0-exp(1.0/epsilon[istate][0])));
    }else if(dim == 2){
        const real x = point[0], y = point[1];
        gradient[0] = (1 + (exp(x/epsilon[istate][0])/epsilon[istate][0])/(1.0-exp(1.0/epsilon[istate][0])))
                    * (y + (exp(y/epsilon[istate][1])-1.0)               /(1.0-exp(1.0/epsilon[istate][1])));
        gradient[1] = (x + (exp(x/epsilon[istate][0])-1.0)               /(1.0-exp(1.0/epsilon[istate][0])))
                    * (1 + (exp(y/epsilon[istate][1])/epsilon[istate][1])/(1.0-exp(1.0/epsilon[istate][1])));
    }else if(dim == 3){
        const real x = point[0], y = point[1], z = point[2];
        gradient[0] = (1 + (exp(x/epsilon[istate][0])/epsilon[istate][0])/(1.0-exp(1.0/epsilon[istate][0])))
                    * (y + (exp(y/epsilon[istate][1])-1.0)               /(1.0-exp(1.0/epsilon[istate][1])))
                    * (z + (exp(z/epsilon[istate][2])-1.0)               /(1.0-exp(1.0/epsilon[istate][2])));
        gradient[1] = (x + (exp(x/epsilon[istate][0])-1.0)               /(1.0-exp(1.0/epsilon[istate][0])))
                    * (1 + (exp(y/epsilon[istate][1])/epsilon[istate][1])/(1.0-exp(1.0/epsilon[istate][1])))
                    * (z + (exp(z/epsilon[istate][2])-1.0)               /(1.0-exp(1.0/epsilon[istate][2])));
        gradient[2] = (x + (exp(x/epsilon[istate][0])-1.0)               /(1.0-exp(1.0/epsilon[istate][0])))
                    * (y + (exp(y/epsilon[istate][1])-1.0)               /(1.0-exp(1.0/epsilon[istate][1])))
                    * (1 + (exp(z/epsilon[istate][2])/epsilon[istate][2])/(1.0-exp(1.0/epsilon[istate][2])));
    }
    return gradient;
}

template <int dim, typename real>
inline dealii::SymmetricTensor<2,dim,real> ManufacturedSolutionSine<dim,real>
::hessian (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    dealii::SymmetricTensor<2,dim,real> hessian;
    // Hard-coded is much more readable than the dimensionally generic one
    const real A = this->amplitudes[istate];
    const dealii::Tensor<1,dim,real> f = this->frequencies[istate];
    if (dim==1) {
        const real fx = f[0]*point[0];
        hessian[0][0] = -A*f[0]*f[0]*sin(fx);
    }
    if (dim==2) {
        const real fx = f[0]*point[0];
        const real fy = f[1]*point[1];
        hessian[0][0] = -A*f[0]*f[0]*sin(fx)*sin(fy);
        hessian[0][1] =  A*f[0]*f[1]*cos(fx)*cos(fy);

        hessian[1][0] =  A*f[1]*f[0]*cos(fx)*cos(fy);
        hessian[1][1] = -A*f[1]*f[1]*sin(fx)*sin(fy);
    }
    if (dim==3) {
        const real fx = f[0]*point[0];
        const real fy = f[1]*point[1];
        const real fz = f[2]*point[2];
        hessian[0][0] = -A*f[0]*f[0]*sin(fx)*sin(fy)*sin(fz);
        hessian[0][1] =  A*f[0]*f[1]*cos(fx)*cos(fy)*sin(fz);
        hessian[0][2] =  A*f[0]*f[2]*cos(fx)*sin(fy)*cos(fz);
        
        hessian[1][0] =  A*f[1]*f[0]*cos(fx)*cos(fy)*sin(fz);
        hessian[1][1] = -A*f[1]*f[1]*sin(fx)*sin(fy)*sin(fz);
        hessian[1][2] =  A*f[1]*f[2]*sin(fx)*cos(fy)*cos(fz);
        
        hessian[2][0] =  A*f[2]*f[0]*cos(fx)*sin(fy)*cos(fz);
        hessian[2][1] =  A*f[2]*f[1]*sin(fx)*cos(fy)*cos(fz);
        hessian[2][2] = -A*f[2]*f[2]*sin(fx)*sin(fy)*sin(fz);
    }
    return hessian;
}

template <int dim, typename real>
inline dealii::SymmetricTensor<2,dim,real> ManufacturedSolutionAdd<dim,real>
::hessian (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    dealii::SymmetricTensor<2,dim,real> hessian;
    const real A = this->amplitudes[istate];
    const dealii::Tensor<1,dim,real> f = this->frequencies[istate];
    if (dim==1) {
        const real fx = f[0]*point[0];
        hessian[0][0] = -A*f[0]*f[0]*sin(fx);
    }
    if (dim==2) {
        const real fx = f[0]*point[0];
        const real fy = f[1]*point[1];
        hessian[0][0] = -A*f[0]*f[0]*sin(fx);
        hessian[0][1] =  0.0;

        hessian[1][0] =  0.0;
        hessian[1][1] = -A*f[1]*f[1]*sin(fy);
    }
    if (dim==3) {
        const real fx = f[0]*point[0];
        const real fy = f[1]*point[1];
        const real fz = f[2]*point[2];
        hessian[0][0] = -A*f[0]*f[0]*sin(fx);
        hessian[0][1] =  0.0;
        hessian[0][2] =  0.0;
        
        hessian[1][0] =  0.0;
        hessian[1][1] = -A*f[1]*f[1]*sin(fy);
        hessian[1][2] =  0.0;
        
        hessian[2][0] =  0.0;
        hessian[2][1] =  0.0;
        hessian[2][2] = -A*f[2]*f[2]*sin(fz);
    }
    return hessian;
}

template <int dim, typename real>
inline dealii::SymmetricTensor<2,dim,real> ManufacturedSolutionCosine<dim,real>
::hessian (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    dealii::SymmetricTensor<2,dim,real> hessian;
    const real A = this->amplitudes[istate];
    const dealii::Tensor<1,dim,real> f = this->frequencies[istate];
    if (dim==1) {
        const real fx = f[0]*point[0];
        hessian[0][0] = -A*f[0]*f[0]*cos(fx);
    }
    if (dim==2) {
        const real fx = f[0]*point[0];
        const real fy = f[1]*point[1];
        hessian[0][0] = -A*f[0]*f[0]*cos(fx)*cos(fy);
        hessian[0][1] =  A*f[0]*f[1]*sin(fx)*sin(fy);

        hessian[1][0] =  A*f[1]*f[0]*sin(fx)*sin(fy);
        hessian[1][1] = -A*f[1]*f[1]*cos(fx)*cos(fy);
    }
    if (dim==3) {
        const real fx = f[0]*point[0];
        const real fy = f[1]*point[1];
        const real fz = f[2]*point[2];
        hessian[0][0] = -A*f[0]*f[0]*cos(fx)*cos(fy)*cos(fz);
        hessian[0][1] =  A*f[0]*f[1]*sin(fx)*sin(fy)*cos(fz);
        hessian[0][2] =  A*f[0]*f[2]*sin(fx)*cos(fy)*sin(fz);
        
        hessian[1][0] =  A*f[1]*f[0]*sin(fx)*sin(fy)*cos(fz);
        hessian[1][1] = -A*f[1]*f[1]*cos(fx)*cos(fy)*cos(fz);
        hessian[1][2] =  A*f[1]*f[2]*cos(fx)*sin(fy)*sin(fz);
        
        hessian[2][0] =  A*f[2]*f[0]*sin(fx)*cos(fy)*sin(fz);
        hessian[2][1] =  A*f[2]*f[1]*cos(fx)*sin(fy)*sin(fz);
        hessian[2][2] = -A*f[2]*f[2]*cos(fx)*cos(fy)*cos(fz);
    }
    return hessian;
}

template <int dim, typename real>
inline dealii::SymmetricTensor<2,dim,real> ManufacturedSolutionExp<dim,real>
::hessian (const dealii::Point<dim,real> &point, const unsigned int  /*istate*/) const
{
    dealii::SymmetricTensor<2,dim,real> hessian;
    if (dim==1) {
        hessian[0][0] = exp(point[0]);
    }
    if (dim==2) {
        hessian[0][0] = exp(point[0]);
        hessian[0][1] = 0.0;

        hessian[1][0] = 0.0;
        hessian[1][1] = exp(point[1]);
    }
    if (dim==3) {
        hessian[0][0] = exp(point[0]);
        hessian[0][1] = 0.0;
        hessian[0][2] = 0.0;
        
        hessian[1][0] = 0.0;
        hessian[1][1] = exp(point[1]);
        hessian[1][2] = 0.0;
        
        hessian[2][0] = 0.0;
        hessian[2][1] = 0.0;
        hessian[2][2] = exp(point[2]);
    }
    return hessian;
}

template <int dim, typename real>
inline dealii::SymmetricTensor<2,dim,real> ManufacturedSolutionEvenPoly<dim,real>
::hessian (const dealii::Point<dim,real> &point, const unsigned int  /*istate*/) const
{
    dealii::SymmetricTensor<2,dim,real> hessian;
    const double poly_max = 7;
    if (dim==1) {
        hessian[0][0] = poly_max*poly_max*pow(point[0] + 0.5, poly_max-2);
    }
    if (dim==2) {
        hessian[0][0] = poly_max*poly_max*pow(point[0] + 0.5, poly_max-2);
        hessian[0][1] = 0.0;

        hessian[1][0] = 0.0;
        hessian[1][1] = poly_max*poly_max*pow(point[1] + 0.5, poly_max-2);
    }
    if (dim==3) {
        hessian[0][0] = poly_max*poly_max*pow(point[0] + 0.5, poly_max-2);
        hessian[0][1] = 0.0;
        hessian[0][2] = 0.0;
        
        hessian[1][0] = 0.0;
        hessian[1][1] = poly_max*poly_max*pow(point[1] + 0.5, poly_max-2);
        hessian[1][2] = 0.0;
        
        hessian[2][0] = 0.0;
        hessian[2][1] = 0.0;
        hessian[2][2] = poly_max*poly_max*pow(point[2] + 0.5, poly_max-2);
    }
    return hessian;
}

template <int dim, typename real>
inline dealii::SymmetricTensor<2,dim,real> ManufacturedSolutionPoly<dim,real>
::hessian (const dealii::Point<dim,real> &point, const unsigned int  /*istate*/) const
{
    dealii::SymmetricTensor<2,dim,real> hessian;
    if (dim==1) {
        const real x = point[0];
        hessian[0][0] = - 2.0 -6*x + 12*x*x - 20*x*x*x + 30*x*x*x*x - 2.500*sin(50*x);
    }
    if (dim==2) {
        real x = point[0];
        hessian[0][0] = - 2.0 -6*x + 12*x*x - 20*x*x*x + 30*x*x*x*x - 2.500*sin(50*x);
        x = point[1];
        hessian[1][1] = - 2.0 -6*x + 12*x*x - 20*x*x*x + 30*x*x*x*x - 2.500*sin(50*x);
    }
    if (dim==3) {
        real x = point[0];
        hessian[0][0] = - 2.0 -6*x + 12*x*x - 20*x*x*x + 30*x*x*x*x - 2.500*sin(50*x);
        x = point[1];
        hessian[1][1] = - 2.0 -6*x + 12*x*x - 20*x*x*x + 30*x*x*x*x - 2.500*sin(50*x);
        x = point[2];
        hessian[2][2] = - 2.0 -6*x + 12*x*x - 20*x*x*x + 30*x*x*x*x - 2.500*sin(50*x);
    }
    return hessian;
}

template <int dim, typename real>
inline dealii::SymmetricTensor<2,dim,real> ManufacturedSolutionAtan<dim,real>
::hessian(const dealii::Point<dim,real> &point, const unsigned int /*istate*/) const
{
    dealii::SymmetricTensor<2,dim,real> hes;

    for(unsigned int k1 = 0; k1 < dim; ++k1){
        // taking the k1^th derivative
        for(unsigned int k2 = 0; k2 < dim; ++k2){
            // taking the k2^th derivative
            real hes_dim = 1;
            for(unsigned int i = 0; i < dim; ++i){
                real x = point[i];
                real val_dim = 0;
                for(unsigned int j = 0; j < n_shocks[i]; ++j){
                    if(i == k1 && i == k2){
                        // taking the second derivative in this dim
                        real coeff = S_j[i][j]*(x-x_j[i][j]);
                        val_dim += -2.0*pow(S_j[i][j],2)*coeff/pow(pow(coeff,2)+1,2);
                    }else if(i == k1 || i == k2){
                        // taking the first derivative in this dim
                        real coeff = S_j[i][j]*(x-x_j[i][j]);
                        val_dim += S_j[i][j]/(pow(coeff,2)+1);
                    }else{
                        // taking the value in this dim
                        val_dim += atan(S_j[i][j]*(x-x_j[i][j]));
                    }
                }
                hes_dim *= val_dim;
            }
            hes[k1][k2] = hes_dim;
        }
    }

    return hes;
}

template <int dim, typename real>
inline dealii::SymmetricTensor<2,dim,real> ManufacturedSolutionBoundaryLayer<dim,real>
::hessian(const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    dealii::SymmetricTensor<2,dim,real> hessian;
    if (dim==1) {
        const real x = point[0];
        hessian[0][0] = (exp(x/epsilon[istate][0])/pow(epsilon[istate][0],2)/(1.0-exp(1.0/epsilon[istate][0])));
    }
    if (dim==2) {
        real x = point[0], y = point[1];
        hessian[0][0] = (exp(x/epsilon[istate][0])/pow(epsilon[istate][0],2)/(1.0-exp(1.0/epsilon[istate][0])))
                      * (y + (exp(y/epsilon[istate][1])-1.0)                /(1.0-exp(1.0/epsilon[istate][1])));
        hessian[0][1] = (1 + (exp(x/epsilon[istate][0])/epsilon[istate][0]) /(1.0-exp(1.0/epsilon[istate][0])))
                      * (1 + (exp(y/epsilon[istate][1])/epsilon[istate][1]) /(1.0-exp(1.0/epsilon[istate][1])));
        
        hessian[1][0] = hessian[0][1];
        hessian[1][1] = (x + (exp(x/epsilon[istate][0])-1.0)                /(1.0-exp(1.0/epsilon[istate][0])))
                      * (exp(y/epsilon[istate][1])/pow(epsilon[istate][1],2)/(1.0-exp(1.0/epsilon[istate][1])));
    }
    if (dim==3) {
        real x = point[0], y = point[1], z = point[2];
        hessian[0][0] = (exp(x/epsilon[istate][0])/pow(epsilon[istate][0],2)/(1.0-exp(1.0/epsilon[istate][0])))
                      * (y + (exp(y/epsilon[istate][1])-1.0)                /(1.0-exp(1.0/epsilon[istate][1])))
                      * (z + (exp(z/epsilon[istate][2])-1.0)                /(1.0-exp(1.0/epsilon[istate][2])));
        hessian[0][1] = (1 + (exp(x/epsilon[istate][0])/epsilon[istate][0]) /(1.0-exp(1.0/epsilon[istate][0])))
                      * (1 + (exp(y/epsilon[istate][1])/epsilon[istate][1]) /(1.0-exp(1.0/epsilon[istate][1])))
                      * (z + (exp(z/epsilon[istate][2])-1.0)                /(1.0-exp(1.0/epsilon[istate][2])));
        hessian[0][2] = (1 + (exp(x/epsilon[istate][0])/epsilon[istate][0]) /(1.0-exp(1.0/epsilon[istate][0])))
                      * (y + (exp(y/epsilon[istate][1])-1.0)                /(1.0-exp(1.0/epsilon[istate][1])))
                      * (1 + (exp(z/epsilon[istate][2])/epsilon[istate][2]) /(1.0-exp(1.0/epsilon[istate][2])));

        hessian[1][0] = hessian[0][1];
        hessian[1][1] = (x + (exp(x/epsilon[istate][0])-1.0)                /(1.0-exp(1.0/epsilon[istate][0])))
                      * (exp(y/epsilon[istate][1])/pow(epsilon[istate][1],2)/(1.0-exp(1.0/epsilon[istate][1])))
                      * (z + (exp(z/epsilon[istate][2])-1.0)                /(1.0-exp(1.0/epsilon[istate][2])));
        hessian[1][2] = (x + (exp(x/epsilon[istate][0])-1.0)                /(1.0-exp(1.0/epsilon[istate][0])))
                      * (1 + (exp(y/epsilon[istate][1])/epsilon[istate][1]) /(1.0-exp(1.0/epsilon[istate][1])))
                      * (1 + (exp(z/epsilon[istate][2])/epsilon[istate][2]) /(1.0-exp(1.0/epsilon[istate][2])));

        hessian[2][0] = hessian[0][2];
        hessian[2][1] = hessian[2][1];
        hessian[2][2] = (x + (exp(x/epsilon[istate][0])-1.0)                /(1.0-exp(1.0/epsilon[istate][0])))
                      * (y + (exp(y/epsilon[istate][1])-1.0)                /(1.0-exp(1.0/epsilon[istate][1])))
                      * (exp(z/epsilon[istate][2])/pow(epsilon[istate][2],2)/(1.0-exp(1.0/epsilon[istate][2])));
    }
    return hessian;
}

template <int dim, typename real>
ManufacturedSolutionFunction<dim,real>
::ManufacturedSolutionFunction (const unsigned int nstate)
    :
    dealii::Function<dim,real>(nstate)
    , nstate(nstate)
    , base_values(nstate)
    , amplitudes(nstate)
    , frequencies(nstate)
{
    const double pi = atan(1)*4.0;
    //const double ee = exp(1);

    for (int s=0; s<(int)nstate; s++) {
        base_values[s] = 1+(s+1.0)/nstate;
        base_values[nstate-1] = 10;
        amplitudes[s] = 0.2*base_values[s]*sin((static_cast<double>(nstate)-s)/nstate);
        for (int d=0; d<dim; d++) {
            //frequencies[s][d] = 2.0 + sin(0.1+s*0.5+d*0.2) *  pi / 2.0;
            frequencies[s][d] = 2.0 + sin(0.1+s*0.5+d*0.2) *  pi / 2.0;
        }
    
    }
}

template <int dim, typename real>
inline dealii::Tensor<1,dim,real> ManufacturedSolutionFunction<dim,real>
::gradient_fd (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    dealii::Tensor<1,dim,real> gradient;
    const double eps=1e-6;
    for (int dim_deri=0; dim_deri<dim; dim_deri++) {
        dealii::Point<dim,real> pert_p = point;
        dealii::Point<dim,real> pert_m = point;
        pert_p[dim_deri] += eps;
        pert_m[dim_deri] -= eps;
        const real value_p = value(pert_p,istate);
        const real value_m = value(pert_m,istate);
        gradient[dim_deri] = (value_p - value_m) / (2*eps);
    }
    return gradient;
}

template <int dim, typename real>
inline dealii::SymmetricTensor<2,dim,real> ManufacturedSolutionFunction<dim,real>
::hessian_fd (const dealii::Point<dim,real> &point, const unsigned int istate) const
{
    dealii::SymmetricTensor<2,dim,real> hessian;
    const double eps=1e-4;
    for (int d1=0; d1<dim; d1++) {
        for (int d2=d1; d2<dim; d2++) {
            dealii::Point<dim,real> pert_p_p = point;
            dealii::Point<dim,real> pert_p_m = point;
            dealii::Point<dim,real> pert_m_p = point;
            dealii::Point<dim,real> pert_m_m = point;

            pert_p_p[d1] += (+eps); pert_p_p[d2] += (+eps);
            pert_p_m[d1] += (+eps); pert_p_m[d2] += (-eps);
            pert_m_p[d1] += (-eps); pert_m_p[d2] += (+eps);
            pert_m_m[d1] += (-eps); pert_m_m[d2] += (-eps);

            const real valpp = value(pert_p_p, istate);
            const real valpm = value(pert_p_m, istate);
            const real valmp = value(pert_m_p, istate);
            const real valmm = value(pert_m_m, istate);

            hessian[d1][d2] = (valpp - valpm - valmp + valmm) / (4*eps*eps);
        }
    }
    return hessian;
}

template <int dim, typename real>
void ManufacturedSolutionFunction<dim,real>
::vector_gradient (
    const dealii::Point<dim,real> &p,
    std::vector<dealii::Tensor<1,dim, real> > &gradients) const
{
    for (unsigned int i = 0; i < nstate; ++i)
        gradients[i] = gradient(p, i);
}


template <int dim, typename real>
inline std::vector<real> ManufacturedSolutionFunction<dim,real>
::stdvector_values (const dealii::Point<dim,real> &point) const
{
    std::vector<real> values(nstate);
    for (unsigned int s=0; s<nstate; s++) { values[s] = value(point, s); }
    return values;
}



template <int dim, typename real>
std::shared_ptr< ManufacturedSolutionFunction<dim,real> > 
ManufacturedSolutionFactory<dim,real>::create_ManufacturedSolution(
    Parameters::AllParameters const *const param, 
    int                                    nstate)
{
    using ManufacturedSolutionEnum = Parameters::ManufacturedSolutionParam::ManufacturedSolutionType;
    ManufacturedSolutionEnum solution_type = param->manufactured_convergence_study_param.manufactured_solution_param.manufactured_solution_type;

    return create_ManufacturedSolution(solution_type, nstate);
}

template <int dim, typename real>
std::shared_ptr< ManufacturedSolutionFunction<dim,real> >
ManufacturedSolutionFactory<dim,real>::create_ManufacturedSolution(
    Parameters::ManufacturedSolutionParam::ManufacturedSolutionType solution_type,
    int                                                                     nstate)
{
    if(solution_type == ManufacturedSolutionEnum::sine_solution){
        return std::make_shared<ManufacturedSolutionSine<dim,real>>(nstate);
    }else if(solution_type == ManufacturedSolutionEnum::cosine_solution){
        return std::make_shared<ManufacturedSolutionCosine<dim,real>>(nstate);
    }else if(solution_type == ManufacturedSolutionEnum::additive_solution){
        return std::make_shared<ManufacturedSolutionAdd<dim,real>>(nstate);
    }else if(solution_type == ManufacturedSolutionEnum::exp_solution){
        return std::make_shared<ManufacturedSolutionExp<dim,real>>(nstate);
    }else if(solution_type == ManufacturedSolutionEnum::poly_solution){
        return std::make_shared<ManufacturedSolutionPoly<dim,real>>(nstate);
    }else if(solution_type == ManufacturedSolutionEnum::even_poly_solution){
        return std::make_shared<ManufacturedSolutionEvenPoly<dim,real>>(nstate);
    }else if(solution_type == ManufacturedSolutionEnum::atan_solution){
        return std::make_shared<ManufacturedSolutionAtan<dim,real>>(nstate);
    }else if(solution_type == ManufacturedSolutionEnum::boundary_layer_solution){
        return std::make_shared<ManufacturedSolutionBoundaryLayer<dim,real>>(nstate);
    }else{
        std::cout << "Invalid Manufactured Solution." << std::endl;
    }

    return nullptr;
}

template class ManufacturedSolutionSine<PHILIP_DIM,double>;
template class ManufacturedSolutionSine<PHILIP_DIM,Sacado::Fad::DFad<double>>;
template class ManufacturedSolutionCosine<PHILIP_DIM,double>;
template class ManufacturedSolutionCosine<PHILIP_DIM,Sacado::Fad::DFad<double>>;
template class ManufacturedSolutionAdd<PHILIP_DIM,double>;
template class ManufacturedSolutionAdd<PHILIP_DIM,Sacado::Fad::DFad<double>>;
template class ManufacturedSolutionExp<PHILIP_DIM,double>;
template class ManufacturedSolutionExp<PHILIP_DIM,Sacado::Fad::DFad<double>>;
template class ManufacturedSolutionPoly<PHILIP_DIM,double>;
template class ManufacturedSolutionPoly<PHILIP_DIM,Sacado::Fad::DFad<double>>;
template class ManufacturedSolutionEvenPoly<PHILIP_DIM,double>;
template class ManufacturedSolutionEvenPoly<PHILIP_DIM,Sacado::Fad::DFad<double>>;
template class ManufacturedSolutionAtan<PHILIP_DIM,double>;
template class ManufacturedSolutionAtan<PHILIP_DIM,Sacado::Fad::DFad<double>>;
template class ManufacturedSolutionBoundaryLayer<PHILIP_DIM,double>;
template class ManufacturedSolutionBoundaryLayer<PHILIP_DIM,Sacado::Fad::DFad<double>>;

template class ManufacturedSolutionFunction<PHILIP_DIM,double>;
template class ManufacturedSolutionFunction<PHILIP_DIM,Sacado::Fad::DFad<double>>;

template class ManufacturedSolutionFactory<PHILIP_DIM,double>;
template class ManufacturedSolutionFactory<PHILIP_DIM,Sacado::Fad::DFad<double>>;
}
