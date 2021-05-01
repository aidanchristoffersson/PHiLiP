
#ifndef __RECONSTRUCT_POLY_H__
#define __RECONSTRUCT_POLY_H__

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/base/polynomial_space.h>

#include <deal.II/grid/tria.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>

#include "physics/manufactured_solution.h"

namespace PHiLiP {

namespace GridRefinement {

/// Enumeration of norms availible to be used in the patchwise reconstruction
enum class NormType{
    H1,
    L2,
    };

// forward declaration of multi-index computation from Dealii
/// Modified implementaiton of polynomial space indexing from Deal.II
/** Computes multi-index of polynomial term for different dimensions. 
  * Based on the Deal.II protected function: 
  * https://www.dealii.org/current/doxygen/deal.II/polynomial__space_8cc_source.html
  */ 
template <int dim>
std::array<unsigned int, dim> compute_index(
    const unsigned int i,
    const unsigned int size);

// funcitons for polynomial reconstruction
template <int dim, int nstate, typename real>
class ReconstructPoly
{

public:
    /// Deleted default constructor
    ReconstructPoly() = delete;

    /// Constructor. Stores required information about the mesh and quadrature rules.
    ReconstructPoly(
        const dealii::hp::DoFHandler<dim>&        dof_handler,           // dof_handler
        const dealii::hp::MappingCollection<dim>& mapping_collection,    // mapping collection
        const dealii::hp::FECollection<dim>&      fe_collection,         // fe collection
        const dealii::hp::QCollection<dim>&       quadrature_collection, // quadrature collection
        const dealii::UpdateFlags&                update_flags);         // update flags for for volume fe

    /// Reinitialze the internal vectors 
    /** These vectors are used to store the obtained derivative
      * values and directions at each mesh element.
      */
    void reinit(const unsigned int n);

    /// Select the Norm to be used in reconstruction
    /** Results in a modified local reconstruction process for the patch of neighbouring elements.
      */ 
    void set_norm_type(const NormType norm_type);

    /// Construct directional derivatives along the chords of the cell
    /** $p+1$ (or rel_order) derivatives are constructed and extracted along the specified directions
      * from the existing cell size. Once all polynomial terms on the surrounding patch are approximated
      * (see reconstruct_norm for description), the derivative components are obtained by evaluating this
      * function along a given direction:
      * 
      * \f[
      *     u_{\bar{\bm{x}}, p}(\bm{x}) = 
      *     \sum_{|\bm{\alpha}| \leq p} {
      *         \frac{\partial^{\bm{\alpha}} u(\bar{\bm{x}})}
      *              {\bm{\alpha}!} 
      *         (\bm{x}-\bar{\bm{x}})^{\bm{\alpha}}
      *         }
      * \f]
      * \f[
      *     D^{p+1}_{\bm{\xi}} u(\bar{\bm{x}}) h^{p+1} 
      *     = u_{\bar{\bm{x}}, p+1}(\bm{x}+h\bm{\xi}) 
      *     - u_{\bar{\bm{x}}, p}(\bm{x}+h\bm{\xi})
      *     = \sum_{i=0}^{p+1}{
      *     \frac{1}
      *          {i! (p+1-i)!} 
      *     \frac{\partial^{p+1} u(\bar{\bm{x}})}
      *          {\partial x^i \partial y^{p+1-i}} 
      *     (x-\bar{x})^i (y-\bar{y})^{p+1-i}}
      * \f]
      * 
      * Ordering is based on the internal dealii numbering. Used in cases where the orientation of the element
      * is not controlled by the refinement procedure (e.g. fixed fraction cases). Gives direct prediction of how
      * error will change with modifying the length of these axes.
      */ 
    void reconstruct_chord_derivative(
        const dealii::LinearAlgebra::distributed::Vector<real>&solution,   ///< Solution approximation to be reconstructed
        const unsigned int                                     rel_order); ///< Relative order of the approximation

    /// Construct the set of largest perpendicular directional derivatives
    /** $p+1$ (or rel_order) derivatives are constructed  and the largest values are extracted. In order to approximate 
      * the "worst case" unit-ball of the error in the high-order case, the method of Dolejsi is used where values are extracted
      * from the maximum direction in the next perpendicular hyperplane to existing directions. In this way, a set of orthogonal directions 
      * is chosen with descending largest derivative orders. In 2D, this process can be written as solving the collection of functions
      * to find maximums based on the patchwise reconstruction of the polynomial (obtained from reconstruction):
      * 
      * \f[
      *     A_1(\bar{\bm{x}}, p) &= \max_{\left\lVert{\bm{\xi}}\right\rVert_2=1}{|D^p_{\bm{\xi}} u(\bar{\bm{x}})|}
      * \f]
      * \f[
      *     \bm{\xi}_1(\bar{\bm{x}}, p) &= \operatornamewithlimits{argmax}_{\left\lVert{\bm{\xi}}\right\rVert_2=1}{|D^p_{\bm{\xi}} u(\bar{\bm{x}})|}
      * \f]
      * \f[
      *     \varphi(\bar{\bm{x}}, p) & \in \left[ 0, 2\pi\right) \quad s.t. \text{ } \bm{\xi}_1 = (cos(\varphi), sin(\varphi))
      * \f]
      * \f[
      *     A_2(\bar{\bm{x}}, p) &= |D^p_{\bm{\xi}_2} u(\bar{\bm{x}})|, \quad \text{where } \bm{\xi}_1 \cdot \bm{\xi}_2 = 0.
      * \f]
      * 
      * In the linear case where $p=2$, these directions and values can be directly extracted from the local hessian reconstruction.
      * Otherewise, for the high-order case, this is performed by sampling an approximately equidistributed set of points to give a "good enough"
      * approximation. In 2D, 180 points are distributed radially on $[0,\pi]$ to give a $1^\circ$ accuracy (second half of angles will be equal 
      * or negative depending on even/odd polynomial order). In 3D, the initial sampling is done using a fibbonaci spiral mapped to the unit sphere 
      * (a fibbonaci sphere, see https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere/26127012#26127012).
      * This provides a good enough apporximation to an even distrubution for this case. $180 \times 180 / 2$ samples are used to give a 
      * similar angular resolution. The second components are extracted from a unit-circle on the perpendicular plane to the largest direction.
      */
    void reconstruct_directional_derivative(
        const dealii::LinearAlgebra::distributed::Vector<real>&solution,   ///< Solution approximation to be reconstructed
        const unsigned int                                     rel_order); ///< Relative order of the approximation

    /// Constructs directional derivates based on the manufactured solution hessian
    /** For $p=2$ only, gets the exact directional derivative components using the spectral decomposition of the hessian:
      * 
      * \f[
      *     H = \left[\begin{matrix} u_{xx} & u_{xy} \\ u_{yx} & u_{yy} \end{matrix}\right] 
      *     = \left[\begin{matrix} \bm{v}_0 & \bm{v}_1 \end{matrix}\right]
      *     \left[\begin{matrix} \lambda_0 & \\ & \lambda_1 \end{matrix}\right]
      *     \left[\begin{matrix} \bm{v}_0^T \\ \bm{v}_1^T \end{matrix}\right]
      * \f]
      * \f[
      *     D^{p=2}_{\bm{\xi}} u(\bar{\bm{x}}) h^{p+1} 
      *     = \bm{\xi}^T H \bm{\xi}
      *     = \sum_{i=0}^{dim} {\lambda_i \left(\bm{\xi}^T \bm{v}_i\right)^2}
      * \f]
      * 
      * Where then $\lambda_i$ (eigenvales) are the directional derivatives and $v_i$ (eigenvectors) are the direction vectors.
      */
    void reconstruct_manufactured_derivative(
        const std::shared_ptr<ManufacturedSolutionFunction<dim,real>>& manufactured_solution, ///< Manufactured solution function
        const unsigned int                                             rel_order);            ///< Relative order of the approximation

private:
    /// Performs polynomial patchwise reconstruction on the current cell in the selected norm
    /** In order to obtain the high-order derivative terms, an enriched polynomial spaced solution $\tilde{u}\in\mathbb{P}^{p+1}$
      * is obtained on the set of neighboring elements, $D(k)$ for the current element $k$. This leads to finding an equality for the
      * inner-product between the original discontinuous solution $u_h$ and the new enriched continuous solution $\tilde{u}$:
      * 
      * \f[
      *    \left< \tilde{u}, \phi \right>_{L \left(D(k)\right)} 
      *    = \left< u_h, \phi \right>_{L \left(D(k)\right)},
      *    \quad \forall \phi \in \mathbb{P}^{p+1}\left(D(k)\right)
      * \f]
      * 
      * Where $L$ is the normed integral space chosen for the reconstruction to take place. This leads to a system of equations
      * for each polynomial shape function on the patch of cells. This is then evaluated on the set of element quadrature 
      * points to a matrix system that can be solved for the coefficients of the enriched solution in the polynomial space:
      * 
      * \f[
      *     \tilde{u}(x,y) = \sum_{i=0}^{N} {a_i \phi_i(x,y)}
      * \f]
      * \f[
      *     \left[\begin{matrix}
      *     <\phi_0,\phi_0>_{N(D(k))} & \cdots & <\phi_0,\phi_N>_{N(D(k))} \\ 
      *     \vdots & \ddots & \vdots \\ 
      *     <\phi_N,\phi_0>_{N(D(k))} & \cdots & <\phi_N,\phi_N>_{N(D(k))}
      *     \end{matrix}\right]
      *     \left[\begin{matrix}
      *     a_0 \\ \vdots \\ a_N \\
      *     \end{matrix}\right]
      *     = \left[\begin{matrix} 
      *     <u_h,\phi_0>_{N(D(k))} \\ \vdots \\ <u_h,\phi_N>_{N(D(k))} 
      *     \end{matrix}\right]
      * \f]
      * 
      * and returned as a vector for further processing of the directional derivatives. For the current class, the polynomial space
      * is selected as the set of non-homogeneous polynomials of maximum order $p+1$. For example, $\phi_i(x,y) = \left[1, x, y, x^2, ...\right]$.
      */ 
    template <typename DoFCellAccessorType>
    dealii::Vector<real> reconstruct_norm(
        const NormType                                          norm_type,
        const DoFCellAccessorType &                             curr_cell,
        const dealii::PolynomialSpace<dim>                      ps,
        const dealii::LinearAlgebra::distributed::Vector<real> &solution);

    /// Performs polynomial patchwise reconstruction on the current cell in the H1 semi-norm
    /** See general form of reconstruct_norm for basic reconstruction problem description. This function
      * is specialized to work based on the selected H1 semi-norm leading to the patchwise norm definition:
      * 
      * \f[
      *     \left<f,g\right>_{H^1(D(k))} 
      *     = \int_{D(k)}{\left[
      *     f(\bm{x})g(\bm{x}) 
      *     + \sum_{i=0}^{dim} {\partial_i f(\bm{x}) \partial_i g(\bm{x})}
      *     \mathrm{d} \bm{x}
      *     \right]}
      * \f]
      * 
      * This requires the polynomial space values and derivatives to be constructed for each pair of shape functions 
      * on the quadrature points and also integrated with the discrete solution approximation values and derivatives.
      */ 
    template <typename DoFCellAccessorType>
    dealii::Vector<real> reconstruct_H1_norm(
        const DoFCellAccessorType &                             curr_cell,
        const dealii::PolynomialSpace<dim>                      ps,
        const dealii::LinearAlgebra::distributed::Vector<real> &solution);

    /// Performs polynomial patchwise reconstruction on the current cell in the L2 norm
    /** See general form of reconstruct_norm for basic reconstruction problem description. This function
      * is specialized to work based on the selected H1 semi-norm leading to the patchwise norm definition:
      * 
      * \f[
      *     \left<f,g\right>_{L^2(D(k))} 
      *     = \int_{D(k)}{\left[
      *     f(\bm{x})g(\bm{x}) 
      *     \mathrm{d} \bm{x}
      *     \right]}
      * \f]
      * 
      * This requires the polynomial space values to be constructed for each pair of shape functions
      * on the quadrature points and also integrated with the discrete solution approximation values.
      */ 
    template <typename DoFCellAccessorType>
    dealii::Vector<real> reconstruct_L2_norm(
        const DoFCellAccessorType &                             curr_cell,
        const dealii::PolynomialSpace<dim>                      ps,
        const dealii::LinearAlgebra::distributed::Vector<real> &solution);

    /// Get the patch of cells surrounding the current cell of DofCellAccessorType
    /** Returns a list of neighbor cells sharing a face (or subface) with the current cell. 
      * Based on dealii::GridTools::get_patch_around_cell and modified to work directly on the dof_handler
      * accesor for use with dealii::hp::DoFHandler instead of needing to be cast back and forth.
      */ 
    template <typename DoFCellAccessorType>
    std::vector<DoFCellAccessorType> get_patch_around_dof_cell(
        const DoFCellAccessorType &cell);

    // member attributes
    const dealii::hp::DoFHandler<dim>&         dof_handler;           ///< Degree of freedom handler for iteration over mesh elements and their nodes
    const dealii::hp::MappingCollection<dim> & mapping_collection;    ///< Collection of mapping rules for reference element conversion
    const dealii::hp::FECollection<dim> &      fe_collection;         ///< Collection of Finite elements to represent discontinuous $hp$ solution space
    const dealii::hp::QCollection<dim> &       quadrature_collection; ///< Collection of quadrature rules used to evaluate volume integrals
    const dealii::UpdateFlags &                update_flags;          ///< Update flags used in obtaining local cell representation

    /// Setting controls the choice of norm used in reconstruction. Set via set_norm_type.
    NormType norm_type;

public:
    /// Derivative values
    /** For each element, array of values indicates the scale of the $p+1$^th (or rel_order) directional
      * derivatives that have been evaluated from the specified chord directions (in reconstruct_chord_derivative)
      * or from the largest orthogonal set of directions (in reconstruct_directional_derivative).
      * These correspond with the numbering of the unit vector directions in derivative_direction. 
      */ 
    std::vector<std::array<real,dim>>                       derivative_value;
    
    /// Derivative directions
    /** For each element, array of unit vectors indicate the direction of the $p+1$^th (or rel_order) directional
      * derivatives that have been evaluated from the specified chord directions (in reconstruct_chord_derivative)
      * or from the largest orthogonal set of directions (in reconstruct_directional_derivative).
      * These correspond with the numbering of the scalar derivative values in derivative_value. 
      */ 
    std::vector<std::array<dealii::Tensor<1,dim,real>,dim>> derivative_direction;

    /// Gets the i^th largest componet of the directional derivative vector as a dealii::Vector
    dealii::Vector<real> get_derivative_value_vector_dealii(
        const unsigned int index);
};

} // namespace GridRefinement

} //namespace PHiLiP

#endif // __RECONSTRUCT_POLY_H__
