.. Label between '.. _' and ':' ; use :ref:`text <label>` for reference
.. _analysis:

*************
Analysis step
*************

This chapter describes the analysis step(s) implemented in LEKF.


Observation representation
==========================

Simulation of observation vector :math:`\mathbf{y}^o` from *true* state vector :math:`\mathbf{x}`:

.. math::
    \mathbf{y}^o ~=~ \mathbf{H}\ \mathbf{x} ~+~ \mathbf{v}
    
with :math:`\mathbf{v}` the *observation represenation* error, assumed to be as sample
drawn from a normal distribution:

.. math::
    \mathbf{v}~\sim~ N(\mathbf{o},\mathbf{R})


Ensemble notation
=================

The ensemble is denoted with:

.. math::
    \mathbf{x}_1, \dots, \mathbf{x}_m

The ensemble mean and covariance are:

.. math::
    \bar{\mathbf{x}} &= \frac{1}{m}\sum_{i=1}^m \mathbf{x}_i \\
    \mathbf{P} &= \frac{1}{m-1}\sum_{i=1}^m (\mathbf{x}_i-\bar{\mathbf{x}})\ (\mathbf{x}_i-\bar{\mathbf{x}})^T

In *ensemble* notation:

.. math::
    \mathbf{X} &=  \left[ \mathbf{x}_1-\bar{\mathbf{x}}, \dots, \mathbf{x}_m-\bar{\mathbf{x}} \right] \\
    \mathbf{P} &= \frac{1}{m-1} \mathbf{X}\ \mathbf{X}^T \\
    \mathbf{Y} &= \mathbf{H}\ \mathbf{X}

In *ensemble square root* notation:

.. math::
    \mathbf{S} &=  \frac{1}{\sqrt{m-1}}\left[ \mathbf{x}_1-\bar{\mathbf{x}}, \dots, \mathbf{x}_m-\bar{\mathbf{x}} \right] \\
    \mathbf{P} &= \mathbf{S}\ \mathbf{S}^T \\
    \mathbf{\Psi} &= \mathbf{H}\ \mathbf{S}


EnKF analysis
=============


Method
------

The EnKF analysis computes a *minimal variance gain* matrix:

.. math::
    \mathbf{K} &= \mathbf{P}\ \mathbf{H}^T\ 
                   \left[\ \mathbf{H}\ \mathbf{P}\ \mathbf{H}^T+\mathbf{R}\ \right]^{-1} \\
               &= \mathbf{S}\ {\mathbf{\Psi}}^T\ 
                   \left[\ \mathbf{\Psi}{\mathbf{\Psi}}^T+\mathbf{R}\ \right]^{-1}

This gain is applied to each member of the ensemble:

.. math::
    \mathbf{x}_i^a ~=~ \mathbf{x}_i 
         ~+~ \mathbf{K}\ \left(\ \mathbf{y}^o\ -\ \mathbf{H}\ \mathbf{x}_i\ +\ \mathbf{v}_i\ \right)

which includes a random sample of the observation representation error:

.. math::
    \mathbf{v}_i~\sim~ N(\mathbf{o},\mathbf{R})


Implementation
--------------

* Compute random observation representation errors for all observations and ensemble members.

* Solve analysis per sub-domain.

  * Collect observations :math:`\mathbf{y}^o`, simulations :math:`\mathbf{\Psi}`,
    and random observation representation errors :math:`\mathbf{v}_i`:
  
    * within own sub-domain
    * observations from other sub-domains that are within a distance of 3.5 :math:`\rho` 
      of the own sub-domain
      
  * Compute the symmetric *innovation covariance* matrix,
    apply covariance localization using Schur product:
  
    .. math::
        \mathbf{\Gamma} ~=~ \left(\mathbf{\Psi}{\mathbf{\Psi}}^T\right)\circ \mathbf{L}+\mathbf{R}
        
  * Factorize using Choleski decomposition to facilitate computation of the gain matrix:
  
    .. math::
        \mathbf{\Gamma} ~=~ \mathbf{L}\ \mathbf{L}^T
  
  * Loop over elements of the state vector
  
    * Solve for these elements the corresponding rows of the gain matrix from:
    
      .. math::
          \mathbf{L}\ \mathbf{L}^T\ \mathbf{K} ~=~ \mathbf{\Psi} {\mathbf{S}}^T
          
    * Analyse the elements of the ensemble states using this gain.


LETKF analysis
==============

The *Local Ensemble Transform Kalman Filter* [Hunt2007]_ provides an analysis that
can be performed per grid cell of the state, based on nearby observations only.
Here we follow the implementation by [Shin2016]_.


Full matrix formulation
-----------------------

The LETKF algortihm uses a different expression for the *gain* matrix.
The same definition is used:

.. math::
    \mathbf{K} &= \mathbf{P}\ \mathbf{H}^T\ 
                   \left[\ \mathbf{H}\ \mathbf{P}\ \mathbf{H}^T+\mathbf{R}\ \right]^{-1} \\

but using:

.. math::
    \left(\ \mathbf{P}\ \mathbf{H}^T\ \mathbf{R}^{-1}\ \mathbf{H} + \mathbf{I}\ \right)\ \mathbf{P}\ \mathbf{H}^T 
      &= \mathbf{P}\ \mathbf{H}^T\ \mathbf{R}^{-1}\ \mathbf{H}\ \mathbf{P}\ \mathbf{H}^T \ +\ \mathbf{P}\ \mathbf{H}^T \\
      &= \mathbf{P}\ \mathbf{H}^T\ \mathbf{R}^{-1}\ \left[\ \mathbf{H}\ \mathbf{P}\ \mathbf{H}^T +\mathbf{R}\ \right]

this is rewritten to:

.. math::
    \mathbf{K}
       &= \left(\ \mathbf{P}\ \mathbf{H}^T\ \mathbf{R}^{-1}\ \mathbf{H} + \mathbf{I}\ \right)^{-1}\ 
          \left(\ \mathbf{P}\ \mathbf{H}^T\ \mathbf{R}^{-1}\ \mathbf{H} + \mathbf{I}\ \right)\ 
          \mathbf{P}\ \mathbf{H}^T\ 
          \left[\ \mathbf{H}\ \mathbf{P}\ \mathbf{H}^T+\mathbf{R}\ \right]^{-1} \\
       &= \left(\ \mathbf{P}\ \mathbf{H}^T\ \mathbf{R}^{-1}\ \mathbf{H} + \mathbf{I}\ \right)^{-1}\ 
          \mathbf{P}\ \mathbf{H}^T\ \mathbf{R}^{-1}\ \left[ \mathbf{H}\ \mathbf{P}\ \mathbf{H}^T + \mathbf{R}\ \right]\ 
          \left[\ \mathbf{H}\ \mathbf{P}\ \mathbf{H}^T+\mathbf{R}\ \right]^{-1} \\
       &= \left(\ \mathbf{P}\ \mathbf{H}^T\ \mathbf{R}^{-1}\ \mathbf{H} + \mathbf{I}\ \right)^{-1}\ 
          \mathbf{P}\ \mathbf{H}^T\ \mathbf{R}^{-1}

This could then be written as:

.. math::
    \mathbf{K} ~=~ \mathbf{P}^a\ \mathbf{H}^T\ \mathbf{R}^{-1}

using a formulation for the *analyzed covariance*:

.. math::
    \mathbf{P}^a ~=~ \left(\ \mathbf{P}\ \mathbf{H}^T\ \mathbf{R}^{-1}\ \mathbf{H} + \mathbf{I}\ \right)^{-1}\ \mathbf{P}  

Note that the analyzed covariance does not depend on the observations, but only one the
forecast covariance :math:`\mathbf{P}` and the observation covariance error :math:`\mathbf{R}`.
The analyzed mean state is:

.. math::
    \mathbf{x}^a ~=~ \mathbf{x} ~+~ \mathbf{P}^a\ \mathbf{H}^T\ \mathbf{R}^{-1}\ \left(\ \mathbf{y}^o\ -\ \overline{\mathbf{H}\ \mathbf{x}_i}\ \right)

Together, :math:`\mathbf{x}^a` and :math:`\mathbf{P}^a` define the analyzed distribution of the state.


Ensemble formulation
--------------------

The equivalent ensemble formulation is:

.. math:: 
    \mathbf{P}^a &= \left[\ \mathbf{Y}^T\ \mathbf{R}^{-1}\ \mathbf{Y}\ +\ (m-1)\ \mathbf{I}\right]^{-1} \\
    \mathbf{w}^a &= \mathbf{P}^a\ \mathbf{Y}^T\ \mathbf{R}^{-1}\ \left(\ \mathbf{y}^o\ -\ \overline{\mathbf{H}\ \mathbf{x}_i}\ \right) \\
    \mathbf{x}^a &= \mathbf{x} ~+~ \mathbf{X}\ \mathbf{w}^a \\
    \mathbf{X}^a &= \mathbf{X}\ \left[(k-1)\mathbf{P}^a\right]^{1/2}
    
Algorithm
---------

The analysis is performed with the following steps.

* Compute ensemble simulations:

  .. math::
      \underset{r}{\mathbf{y}_i} ~=~ \underset{r\times n}{\mathbf{H}}\ \underset{n}{\mathbf{x}_i}

* Compute ensemble means:

  .. math::
      \mathbf{\bar{x}} &= \overline{\mathbf{x}_i} \\
      \mathbf{\bar{y}} &= \overline{\mathbf{y}_i}
      
* Fill ensemble matrices:

  .. math::
      \underset{n\times m}{\mathbf{X}} &= \left[\ \dots,\mathbf{x}_i-\mathbf{\bar{x}},\dots\ \right] \\
      \underset{r\times m}{\mathbf{Y}} &= \left[\ \dots,\mathbf{y}_i-\mathbf{\bar{y}},\dots\ \right]

* Compute:

  .. math::
      \underset{m\times r}{\mathbf{C}} ~=~ \mathbf{Y}^T\ \mathbf{R}^{-1}
      
  Eventually decrease influcence of observations further away using a localization matrix:
  
  .. math::
      \underset{m\times r}{\mathbf{C}} ~:=~ \underset{m\times r}{\mathbf{C}}\circ\underset{r\times r}{\mathbf{L}}

* Compute eigenvalues and vectors from symmetric positive definite matrix:

  .. math:: 
      \underset{m\times m}{\left[\ \mathbf{C}\ \mathbf{Y}\ +\ (m-1)\ \mathbf{I}\right]}\ \underset{m\times m}{\mathbf{Q}}
         ~=~ \underset{m\times m}{\mathbf{Q}}\ \underset{m\times m}{\mathbf{\Lambda}}
         
  Use that :math:`\mathbf{Q}` is orhonormal and therefore:
  
  .. math::
      \mathbf{Q}^{T}\ \mathbf{Q} &= \mathbf{I} \\
      \mathbf{Q}^{-1} &= \mathbf{Q}^T

  to obtain:
     
  .. math:: 
      \left[\ \mathbf{C}\ \mathbf{Y}\ +\ (m-1)\ \mathbf{I}\right]
         ~=~ \mathbf{Q}\ \mathbf{\Lambda}\ \mathbf{Q}^T

* Compute the analyzed covariance in ensemble space:
   
  .. math:: 
    \underset{m\times m}{\mathbf{P}^a} ~=~ \mathbf{Q}\ \mathbf{\Lambda}^{-1}\ \mathbf{Q}^T
  
  
* Compute the symmetric factorization:
   
  .. math::
      (m-1)\ \mathbf{P}^a ~=~ \mathbf{W}^a\ {\mathbf{W}^a}^T
      
  Use that:
  
  .. math::
      (m-1)\ \mathbf{P}^a &= (m-1)\ \mathbf{Q}\ \mathbf{\Lambda}^{-1}\ \mathbf{Q}^T \\
          &= (m-1)\ \mathbf{Q}\ \mathbf{\Lambda}^{-1/2}\ \mathbf{Q}^T\ \mathbf{Q}\ \mathbf{\Lambda}^{-1/2}\ \mathbf{Q}^T
        
  to obtain:
  
  .. math::
      \underset{m\times m}{\mathbf{W}^a} ~=~ \sqrt{m-1}\ \mathbf{Q}\ \mathbf{\Lambda}^{-1/2}\ \mathbf{Q}^T

* Compute analysis weights:

  .. math::
      \underset{m}{\mathbf{\bar{w}}^a} ~=~ \underset{m\times m}{\mathbf{P}^a}\ \underset{m\times r}{\mathbf{C}}
      \ \underset{r}{\left(\ \mathbf{y}^o\ -\ \mathbf{\bar{y}}\ \right)}

* Compute analysis mean:

  .. math::
      \mathbf{\bar{x}}^a ~=~ \mathbf{\bar{x}} ~+~ \mathbf{X}\ \mathbf{\bar{w}}^a

* Compute analysis ensemble using the columns :math:`\mathbf{w}^a_i` of :math:`\mathbf{W}`:

  .. math::
      \mathbf{x}^a_i ~=~ \mathbf{\bar{x}}^a ~+~ \mathbf{X}\ \mathbf{w}^a_i

  or equivalent:
  
  .. math::
      \mathbf{x}^a_i ~=~ \mathbf{\bar{x}} ~+~ \mathbf{X}\ \left(\mathbf{\bar{w}}^a+\mathbf{w}^a_i\right)


Local analysis
--------------

In the LETKF algorithm, the above analysis is applied per grid cell. 
The algorithm becomes:

* Compute in each domain simulated observations for all ensemble members.
* Collect per domain also the observations from neighbouring domains that are within :math:`3.5\rho` distance
* Loop over blocks in the state: concentrations, aerosol water, noise factors

  * For each block, loop over grid cells

    * Select observations :math:`\tilde{\mathbf{y}}^o` and simulations
      :math:`\tilde{\mathbf{y}}_i` that are:
      
      * within range  :math:`3.5\rho` of the grid cell;
      * specified in the settings as to be used to analyze this block;

    * Compute analysis weights, use localization correlation with decay with distance.    
    * Apply the analysis with the ensemble elements for the selected grid cell 
      and the selected (simulated) observations




