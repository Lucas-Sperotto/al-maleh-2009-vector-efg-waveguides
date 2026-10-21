"""Synthetic and archival contract tests for preliminary cross-method vector fields."""
import importlib.util
import math
from pathlib import Path
import unittest
import numpy as np

SOURCE=Path(__file__).resolve().parents[1]/"scripts"/"m6f_centroid_overlap.py"
spec=importlib.util.spec_from_file_location("m6f",SOURCE)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class M6fFieldTests(unittest.TestCase):
    def test_sign_and_scale_do_not_change_correlation(self):
        a=np.array([1.,-2.,3.,-4.])
        w=np.array([.25,.25,.25,.25])
        self.assertAlmostEqual(mod.correlation(a,-2*a,w),1,places=12)
        self.assertAlmostEqual(mod.correlation(a,2*a,w),1,places=12)

    def test_zero_mode_fails_explicitly(self):
        a=np.zeros(4)
        with self.assertRaises(ValueError):
            mod.correlation(a,a,np.ones(4))

    def test_orthogonal_vectors(self):
        a=np.array([1.,0.,0.,0.])
        b=np.array([0.,1.,0.,0.])
        self.assertAlmostEqual(mod.correlation(a,b,np.ones(4)),0,places=12)

    def test_doublet_basis_rotation_invariance(self):
        rng=np.random.default_rng(7)
        A=rng.normal(size=(50,2))
        B=A @ np.array([[.8,.6],[-.6,.8]])
        weights=np.ones(50)/25
        cos=mod.principal_correlations(A,B,weights)
        np.testing.assert_allclose(cos,[1.,1.],atol=1e-12)
        B2=B @ np.diag([5.,-2.])
        np.testing.assert_allclose(
            mod.principal_correlations(A,B2,weights),cos,atol=1e-12)

    def test_doublet_nontrivial_angle(self):
        A=np.array([[1.,0.],[0.,1.],[0.,0.],[0.,0.]])
        B=np.array([[1.,0.],[0.,math.cos(math.pi/6)],
                    [0.,0.],[0.,math.sin(math.pi/6)]])
        c=mod.principal_correlations(A,B,np.ones(4))
        np.testing.assert_allclose(c,[1.,math.cos(math.pi/6)],atol=1e-12)

    def test_rank_deficient_doublet_is_rejected(self):
        A=np.array([[1.,1.],[0.,0.],[0.,0.],[0.,0.]])
        with self.assertRaises(ValueError):
            mod.principal_correlations(A,A,np.ones(4))

    def test_reference_vector_H_is_divergence_free_analytically(self):
        x=np.array([0.125,0.625])
        y=np.array([0.25,0.75])
        H=mod.analytic_H(x,y,1,1)
        kx=math.pi/2
        ky=math.pi
        expected_x=-ky*np.sin(kx*x)*np.cos(ky*y)
        expected_y=kx*np.cos(kx*x)*np.sin(ky*y)
        np.testing.assert_allclose(H[0::2],expected_x)
        np.testing.assert_allclose(H[1::2],expected_y)

if __name__=="__main__":
    unittest.main()
