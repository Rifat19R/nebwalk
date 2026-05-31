"""
Tests for nebwalk.interpolate.idpp_interpolate.

Run with:  pytest tests/test_idpp.py -v
"""

import numpy as np
import pytest
from ase import Atoms
from nebwalk.interpolate import (
    idpp_interpolate,
    linear_interpolate,
    _pairwise,
    _idpp_obj_and_grad,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _h2(x_shift=0.0):
    """H2 at 0.74 Å bond length, optionally shifted along x."""
    return Atoms('H2', positions=[[x_shift, 0.0, 0.0],
                                  [x_shift + 0.74, 0.0, 0.0]])


def _ethane_staggered():
    """Approximate staggered ethane (phi ~ 60°)."""
    pos = np.array([
        [ 0.000,  0.000,  0.000],   # C1
        [ 1.540,  0.000,  0.000],   # C2
        [-0.390,  1.027,  0.000],   # H1
        [-0.390, -0.513,  0.889],   # H2
        [-0.390, -0.513, -0.889],   # H3
        [ 1.930,  1.027,  0.000],   # H4
        [ 1.930, -0.513,  0.889],   # H5
        [ 1.930, -0.513, -0.889],   # H6
    ])
    return Atoms('C2H6', positions=pos)


def _ethane_eclipsed():
    """Approximate eclipsed ethane (phi ~ 0°)."""
    pos = np.array([
        [ 0.000,  0.000,  0.000],
        [ 1.540,  0.000,  0.000],
        [-0.390,  1.027,  0.000],
        [-0.390, -0.513,  0.889],
        [-0.390, -0.513, -0.889],
        [ 1.930,  0.513,  0.889],   # rotated 60°
        [ 1.930,  0.513, -0.889],
        [ 1.930, -1.027,  0.000],
    ])
    return Atoms('C2H6', positions=pos)


# ---------------------------------------------------------------------------
# Output shape and type
# ---------------------------------------------------------------------------

class TestOutputShape:
    def test_length(self):
        imgs = idpp_interpolate(_h2(), _h2(2.0), n_images=5)
        assert len(imgs) == 7

    def test_endpoints_preserved(self):
        start, end = _h2(), _h2(2.0)
        imgs = idpp_interpolate(start, end, n_images=3)
        np.testing.assert_allclose(imgs[0].positions, start.positions, atol=1e-10)
        np.testing.assert_allclose(imgs[-1].positions, end.positions, atol=1e-10)

    def test_n_atoms_unchanged(self):
        imgs = idpp_interpolate(_h2(), _h2(2.0), n_images=4)
        for img in imgs:
            assert len(img) == 2

    def test_atom_symbols_preserved(self):
        imgs = idpp_interpolate(_h2(), _h2(2.0), n_images=3)
        for img in imgs:
            assert list(img.symbols) == ['H', 'H']

    def test_no_calculator_on_images(self):
        imgs = idpp_interpolate(_h2(), _h2(2.0), n_images=3)
        for img in imgs:
            assert img.calc is None

    def test_single_image(self):
        start, end = _h2(), _h2(2.0)
        imgs = idpp_interpolate(start, end, n_images=1)
        assert len(imgs) == 3


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

class TestValidation:
    def test_species_mismatch_raises(self):
        start = Atoms('H2', positions=np.zeros((2, 3)))
        end   = Atoms('He2', positions=np.ones((2, 3)))
        with pytest.raises(ValueError, match="species"):
            idpp_interpolate(start, end, n_images=2)

    def test_length_mismatch_raises(self):
        start = Atoms('H2', positions=np.zeros((2, 3)))
        end   = Atoms('H3', positions=np.zeros((3, 3)))
        with pytest.raises(ValueError, match="number of atoms"):
            idpp_interpolate(start, end, n_images=2)

    def test_n_images_zero_raises(self):
        with pytest.raises(ValueError):
            idpp_interpolate(_h2(), _h2(2.0), n_images=0)

    def test_pbc_raises(self):
        from ase.build import bulk
        start = bulk('Al', 'fcc', a=4.05, cubic=True)
        end   = start.copy()
        end.positions[0] += [0.5, 0.5, 0.0]
        with pytest.raises(NotImplementedError, match="periodic"):
            idpp_interpolate(start, end, n_images=3)


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------

class TestPhysics:
    def test_bond_length_preserved_for_pure_translation(self):
        """
        For a pure rigid translation, both endpoints have the same internal
        geometry. IDPP target distances are identical at every image, so the
        optimised positions must preserve the bond length.
        """
        start = _h2(0.0)
        end   = _h2(3.0)   # same geometry, shifted 3 Å
        imgs  = idpp_interpolate(start, end, n_images=7)
        for img in imgs:
            _, dist = _pairwise(img.positions)
            bond = dist[0, 1]
            assert abs(bond - 0.74) < 0.01, \
                f"Bond length deviated to {bond:.4f} Å (expected 0.74 Å)"

    def test_different_from_linear_for_torsion(self):
        """
        For torsional motion, IDPP must produce different intermediate positions
        from linear interpolation. If they are identical, IDPP is doing nothing.
        """
        start = _ethane_staggered()
        end   = _ethane_eclipsed()
        n = 5
        imgs_lin  = linear_interpolate(start, end, n_images=n)
        imgs_idpp = idpp_interpolate(start, end, n_images=n)
        any_different = any(
            np.max(np.abs(imgs_lin[k].positions - imgs_idpp[k].positions)) > 1e-4
            for k in range(1, n + 1)
        )
        assert any_different, "IDPP produced identical positions to linear interpolation"

    def test_lower_idpp_objective_than_linear(self):
        """
        The IDPP path must have a strictly lower IDPP objective value than the
        linear path. This is the fundamental guarantee of the method.
        """
        start = _ethane_staggered()
        end   = _ethane_eclipsed()
        n = 7
        n_total = n + 2

        _, d_start = _pairwise(start.positions)
        _, d_end   = _pairwise(end.positions)

        def total_obj(imgs):
            total = 0.0
            for k in range(1, n_total - 1):
                alpha    = k / (n_total - 1)
                d_target = d_start + alpha * (d_end - d_start)
                with np.errstate(divide='ignore', invalid='ignore'):
                    weights = np.where(d_target > 1e-10, 1.0 / d_target ** 4, 0.0)
                np.fill_diagonal(weights, 0.0)
                obj, _ = _idpp_obj_and_grad(
                    imgs[k].positions.ravel(), len(start), d_target, weights
                )
                total += obj
            return total

        imgs_lin  = linear_interpolate(start, end, n_images=n)
        imgs_idpp = idpp_interpolate(start, end, n_images=n)

        obj_lin  = total_obj(imgs_lin)
        obj_idpp = total_obj(imgs_idpp)
        assert obj_idpp < obj_lin, \
            f"IDPP objective ({obj_idpp:.4f}) >= linear ({obj_lin:.4f})"


# ---------------------------------------------------------------------------
# Gradient correctness (finite difference)
# ---------------------------------------------------------------------------

class TestGradient:
    def test_analytical_matches_finite_difference(self):
        """
        The analytical gradient of _idpp_obj_and_grad must match the central
        finite difference to rtol=1e-4. This is the primary numerical safeguard.
        """
        rng = np.random.default_rng(42)
        n_atoms = 4
        positions = rng.standard_normal((n_atoms, 3))

        d_target = np.abs(rng.standard_normal((n_atoms, n_atoms))) + 0.5
        d_target = 0.5 * (d_target + d_target.T)
        np.fill_diagonal(d_target, 0.0)

        with np.errstate(divide='ignore'):
            weights = np.where(d_target > 1e-10, 1.0 / d_target ** 4, 0.0)
        np.fill_diagonal(weights, 0.0)

        x0 = positions.ravel()
        _, grad_analytical = _idpp_obj_and_grad(x0, n_atoms, d_target, weights)

        eps = 1e-5
        grad_fd = np.zeros_like(x0)
        for i in range(len(x0)):
            xp, xm = x0.copy(), x0.copy()
            xp[i] += eps
            xm[i] -= eps
            fp, _ = _idpp_obj_and_grad(xp, n_atoms, d_target, weights)
            fm, _ = _idpp_obj_and_grad(xm, n_atoms, d_target, weights)
            grad_fd[i] = (fp - fm) / (2 * eps)

        np.testing.assert_allclose(
            grad_analytical, grad_fd, rtol=1e-4, atol=1e-6,
            err_msg="Analytical gradient does not match finite differences"
        )
