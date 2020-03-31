#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


class FiniteElement:
    '''
    A class used to perform finite element analysis

    Attributes
    ----------
    problem_type : str
        the identifier of the finite element problem type
    '''

    def __init__(self, problem_type='Elasticity'):
        '''
        Parameters
        ----------
        problem_type : str
        finite element problem type identifier
        '''

        self.problem_type = problem_type

    def get_Ke(self, Ae=None, Ee=None, le=None, EIe=None):
        '''Return element stiffness matrix
        Parameters
        ----------
        Ae : double
            Area of element,  for elasticity
        Ee : double
            Young's modulus of element,  for elasticity
        le : double
            Length of element
        EIe : double
            Product of Young's modulus and moment of inertia,  for beam theory
        '''

        if self.problem_type == 'Elasticity':
            return (Ae * Ee / le) * np.array([[1., - 1.], [-1., 1.]])
        elif self.problem_type == 'EB_beam':
            ee = EIe/(le**3)*np.array([12,6*le,-12,6*le,6*le,4*le**2,-6*le,2*le**2,-12,-6*le,12,-6*le,6*le,2*le**2,-6*le,4*le**2])
            return ee
        else:
            raise Exception('Invalid choice of element type')

    def get_fe_omega(self, le=None, b1=None, b2=None, fe=None):
        '''Return force vector due to distributed loading
        Parameters
        ----------
        le : double
            Length of element
        b1 : double
            Body force at node 1, for elasticity
        b2 : double
            Body force at node 2, for elasticity
        fe : double
           Average distributed load on element
        '''
        if self.problem_type == 'Elasticity':
            return (le / 6.) * np.array([[2., 1.], [1., 2.]]).dot(np.array([[b1], [b2]]))
        elif self.problem_type == 'EB_beam':
            ff = (fe*le/2.0)*np.array([1.0,le/6.0,1.0,-le/6.0])
            return ff
        else:
            raise Exception('Invalid choice of element type')

    def get_dynamic_fe_omega(self, le=None, b1=None, b2=None, fe=None,node= None):
        if self.problem_type == 'EB_beam':
            ff = (fe[node]*le/2.0)*np.array([1.0,le/6.0,1.0,-le/6.0])
            return ff
        else:
            raise Exception('Invalid choice of element type')

    def get_gauss_quadrature_index(self, f = None, a = None, b =None):
        return 0.5*(a+b)+0.5*(b-a)*f
                   
    def get_dof_index(self, e):
        '''Return the global dof associated with an element
        Parameters
        ----------
        e : int
            Element ID
        '''
        if self.problem_type == 'Elasticity':
            return [e, e + 1]
        elif self.problem_type == 'EB_beam':
            exy = np.zeros((2,16),dtype=np.int)
            ey = np.zeros(16)
            for ii in range(4):
                for jj in range(4):
                    exy[0,ii*4+jj] = e*2 +ii;
            for ii in range(4):
                for jj in range(4):
                    exy[1,ii*4+jj] = e*2 +jj;
            return exy
        else:
            raise Exception('Invalid choice of element type')

    def get_strain_e(self, le, de):
        '''Return the strain
        Parameters
        ----------
        le : double
            Length of element
        de : 1D numpy array
            Array of nodal displacements for element
        '''
        if self.problem_type == 'Elasticity':
            return (1. / le) * np.array([-1., 1.]).dot(de)
        else:
            raise Exception('Invalid choice of element type')
