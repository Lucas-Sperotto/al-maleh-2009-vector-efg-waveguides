"""Scientific invariants for reusable independent-validation postprocessing."""
import importlib.util
import csv
import copy
import gzip
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy import io

MODULE_PATH = Path(__file__).resolve().parents[1] / 'scripts' / 'independent_validation.py'
SPEC = importlib.util.spec_from_file_location('independent_validation', MODULE_PATH)
v = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(v)


class IndependentValidationTests(unittest.TestCase):
    def setUp(self):
        x, y = np.meshgrid((np.arange(40) + .5) * 2 / 40, (np.arange(20) + .5) / 20)
        self.points = np.column_stack((x.ravel(), y.ravel()))
        self.weights = np.full(len(self.points), 2 / len(self.points))

    def test_analytic_eigenvalue_degeneracy_and_units(self):
        self.assertAlmostEqual(v.eigenvalue(2, 2), v.eigenvalue(4, 1))
        self.assertAlmostEqual(v.eigenvalue(1, 1), 1.25 * np.pi ** 2)
        self.assertAlmostEqual(np.sqrt(v.eigenvalue(2, 1)) ** 2, 2 * np.pi ** 2)

    def test_rotated_degenerate_subspace_is_preserved(self):
        fields = v.analytic_fields(self.points, ((2, 2), (4, 1)))
        theta = .713
        rotation = np.array([[np.cos(theta), -np.sin(theta)],
                             [np.sin(theta), np.cos(theta)]])
        metrics = v.principal_metrics(fields, fields @ rotation, self.weights)
        self.assertGreater(metrics['minimum_correlation'], 1 - 1e-14)
        self.assertLess(metrics['projector_distance_frobenius'], 1e-7)
        # Individual vector correlation changes under this valid basis rotation.
        original = v.normalize_fields(fields, self.weights)
        rotated = v.normalize_fields(fields @ rotation, self.weights)
        self.assertLess(abs(np.dot(self.weights * original[:, 0], rotated[:, 0])), .8)

    def test_assignment_uses_fields_when_spectral_nearest_is_wrong(self):
        fields = v.analytic_fields(self.points)
        values = np.array([v.eigenvalue(*mode) for mode in v.MODES])
        # Reverse the simple TM31/TM12 eigenvalue locations (both remain in window).
        values[[2, 3]] = values[[3, 2]]
        fields[:, 4:6] = fields[:, 4:6] @ np.array([[.8, -.6], [.6, .8]])
        assignment, _, _ = v.assign_modes(values, fields, self.points, self.weights)
        self.assertEqual(assignment[(3, 1)], [2])
        self.assertEqual(assignment[(1, 2)], [3])
        self.assertEqual(set(assignment[(2, 2)]), {4, 5})
        self.assertEqual(assignment[(2, 2)], assignment[(4, 1)])

    def test_rank_deficient_subspace_is_rejected(self):
        fields = v.analytic_fields(self.points, [(1, 1)])
        with self.assertRaises(ValueError):
            v.principal_metrics(np.hstack((fields, fields)), np.hstack((fields, fields)), self.weights)

    def test_barycentric_interpolation_uses_triangle_basis(self):
        nodes = np.array([[0., 0.], [2., 0.], [2., 1.], [0., 1.]])
        elements = np.array([[0, 1, 2], [0, 2, 3]])
        interpolation = v.p1_interpolation(nodes, elements, np.arange(4), self.points)
        nodal_values = 2 + 3 * nodes[:, 0] - 4 * nodes[:, 1]
        np.testing.assert_allclose(interpolation @ nodal_values,
                                   2 + 3 * self.points[:, 0] - 4 * self.points[:, 1], atol=1e-14)

    def test_eliminated_dirichlet_nodes_are_exact_zero(self):
        nodes = np.array([[0., 0.], [2., 0.], [1., 1.]])
        result = v.p1_interpolation(nodes, np.array([[0, 1, 2]]), np.array([2]), np.array([[1., 0.], [1., .5]]))
        np.testing.assert_allclose(result @ np.array([2.]), [0, 1])

    def test_nonfinite_artifact_is_rejected(self):
        for invalid in (np.nan, np.inf, -np.inf):
            with self.assertRaises(ValueError):
                v.finite('test', np.array([invalid]))

    def test_residual_tracks_eigenvector_column(self):
        k, m = np.diag([0., 2., 5.]), np.eye(3)
        exact = v.pair_residuals(k, m, np.array([0., 2., 5.]), np.eye(3))
        swapped = v.pair_residuals(k, m, np.array([0., 2., 5.]), np.eye(3)[:, [0, 2, 1]])
        self.assertEqual(exact.max(), 0)
        self.assertGreater(swapped.max(), .1)

    def test_exact_same_matrices_lapack_retains_kernel_and_doublet(self):
        k, m = np.diag([0., 0., 2., 2., 5.]), np.diag([1., 2., 1., 1., 1.])
        vectors = np.diag(1 / np.sqrt(np.diag(m)))
        problem = dict(k=k, m=m, values=np.array([0., 0., 2., 2., 5.]), vectors=vectors,
                       residuals=np.zeros(5), cutoff=5e-9)
        rows, clusters, summary = v.solver_comparison(problem, 'synthetic')
        self.assertEqual(len(rows), 5)
        self.assertEqual(clusters[0]['classification'], 'numerical_kernel')
        self.assertEqual(clusters[0]['dimension'], 2)
        self.assertEqual([r['dimension'] for r in clusters], [2, 2, 1])
        self.assertEqual(summary['backend'], 'LAPACKE_dsygv')

    def test_roundtrip_complete_spectrum_and_m_normalization(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            k, m = np.diag([0., 2., 5.]), np.diag([2., 1., .5])
            vectors = np.diag(1 / np.sqrt(np.diag(m)))
            values = np.array([0., 2., 10.])
            for name, matrix in [('stiffness_matrix', k), ('mass_matrix', m), ('eigenvectors', vectors)]:
                with gzip.open(path / (name + '.mtx.gz'), 'wb') as stream:
                    io.mmwrite(stream, matrix, precision=17)
            v.write_csv(path / 'eigenvalues.csv', [dict(index=i, **{'lambda': x}) for i, x in enumerate(values)])
            (path / 'run_metadata.json').write_text(json.dumps({'geometry': {'a': 2, 'b': 1}, 'family': 'TM', 'boundary_condition': 'homogeneous Dirichlet Ez=0'}))
            problem = v.audit_problem(path, 'fem')
            self.assertEqual(problem['audit']['complete_vectors'], 3)
            self.assertEqual(problem['audit']['numerical_kernel_dimension'], 1)
            self.assertLess(problem['audit']['m_orthogonality_max_error'], 1e-14)
            np.testing.assert_array_equal(problem['values'], values)
            metadata = json.loads((path / 'run_metadata.json').read_text())
            metadata['geometry'] = {}
            (path / 'run_metadata.json').write_text(json.dumps(metadata))
            with self.assertRaisesRegex(ValueError, 'explicitly recorded'):
                v.audit_problem(path, 'fem')

    def test_roundtrip_rejects_index_column_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            v.write_csv(path / 'eigenvalues.csv', [dict(index=1, **{'lambda': 1})])
            with self.assertRaisesRegex(ValueError, 'column convention'):
                v.audit_problem(path, 'fem')

    def test_mutual_field_overlap_is_required_for_cross_method_success(self):
        exact = v.analytic_fields(self.points)
        contaminant = v.analytic_fields(self.points, [(5, 1)])[:, 0]
        problems = []
        for sign in (1, -1):
            fields = exact.copy()
            fields[:, 0] += sign * .75 * contaminant
            values = np.array([v.eigenvalue(*mode) for mode in v.MODES])
            assignment, normalized, overlaps = v.assign_modes(values, fields, self.points, self.weights)
            problems.append(dict(values=values, assignments=assignment, fields=fields,
                                 normalized=normalized, overlaps=overlaps, residuals=np.zeros(6),
                                 diagnostics=[], audit={'dirichlet_boundary_exact_zero': True}))
        problems[0]['diagnostics'] = [dict(divergence_ratio=0, normal_ratio=0) for _ in v.MODES]
        rows = v.compare(*problems, self.points, self.weights)
        self.assertEqual(rows[0]['efg_status'], 'supported_candidate')
        self.assertEqual(rows[0]['fem_status'], 'supported_candidate')
        self.assertLess(rows[0]['field_correlation'], v.THRESHOLDS['minimum_overlap'])
        self.assertEqual(rows[0]['status'], 'physical_validation_failed')

    def test_report_audit_tampering_is_detected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'report_data.json'
            expected = dict(audits=[dict(dofs=126, status='pass')], environment={'python': 'old'})
            path.write_text(json.dumps(expected))
            v.verify_report(path, dict(audits=[dict(dofs=126, status='pass')], environment={'python': 'new'}))
            with self.assertRaises(ValueError):
                v.verify_report(path, dict(audits=[dict(dofs=510, status='pass')]))

    def solver_archive_fixture(self, path):
        # Valid archived measurements may differ from a new BLAS execution.
        problem = dict(k=np.diag([0., 0., 2., 2., 5.]), m=np.eye(5),
                       values=np.array([0., 0., 2., 2., 5.]), vectors=np.eye(5),
                       residuals=np.zeros(5), cutoff=5e-9)
        rows, clusters, summary = v.solver_comparison(problem, 'synthetic')
        archived_rows, archived_clusters = copy.deepcopy(rows), copy.deepcopy(clusters)
        archived_summary = copy.deepcopy(summary)
        for row in archived_rows:
            row['lapack_lambda'] = row['eigen_lambda'] + .5 * row['allowed_absolute_difference']
            row['absolute_difference'] = abs(row['eigen_lambda'] - row['lapack_lambda'])
            row['relative_difference'] = row['absolute_difference'] / max(
                abs(row['eigen_lambda']), abs(row['lapack_lambda']), problem['cutoff'])
            row['allowed_absolute_difference'] = (summary['eigenvalue_absolute_tolerance'] +
                summary['eigenvalue_relative_tolerance'] * max(abs(row['eigen_lambda']), abs(row['lapack_lambda'])))
            row['lapack_relative_residual'] = 5e-8
        for cluster in archived_clusters:
            cluster['minimum_mass_correlation'] = 1 - 4e-7
            cluster['maximum_principal_angle_degrees'] = float(np.degrees(np.arccos(1 - 4e-7)))
        archived_summary.update(
            library='libdifferentlapack.dylib',
            maximum_eigenvalue_absolute_difference=max(row['absolute_difference'] for row in archived_rows),
            maximum_relative_residual=5e-8, m_orthogonality_max_error=2e-6,
            minimum_subspace_mass_correlation=1 - 4e-7)
        v.write_csv(path / 'solver_comparison.csv', archived_rows)
        v.write_csv(path / 'solver_subspaces.csv', archived_clusters)
        (path / 'report_data.json').write_text(json.dumps({'solver_comparison': {'synthetic': archived_summary}}))
        return rows, clusters, {'synthetic': summary}, {'synthetic': problem['cutoff']}

    def test_solver_archive_valid_backend_noise_and_library_difference_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            expected = self.solver_archive_fixture(path)
            archived = v.verify_solver_publication(path, *expected)
            v.verify_report(path / 'report_data.json', {'solver_comparison': expected[2]}, solver_archive=archived)
            self.assertEqual(archived['synthetic']['library'], 'libdifferentlapack.dylib')

    def test_solver_archive_eigenvalue_beyond_original_tolerance_fails(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            expected = self.solver_archive_fixture(path)
            rows = v.read_csv(path / 'solver_comparison.csv')
            rows[0]['lapack_lambda'] = 2 * float(rows[0]['allowed_absolute_difference'])
            v.write_csv(path / 'solver_comparison.csv', rows)
            with self.assertRaisesRegex(ValueError, 'eigenvalue exceeds declared tolerance'):
                v.verify_solver_publication(path, *expected)

    def test_solver_archive_declared_tolerance_cannot_be_increased(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            expected = self.solver_archive_fixture(path)
            rows = v.read_csv(path / 'solver_comparison.csv')
            rows[-1]['allowed_absolute_difference'] = 2 * float(rows[-1]['allowed_absolute_difference'])
            v.write_csv(path / 'solver_comparison.csv', rows)
            with self.assertRaisesRegex(ValueError, 'allowed_absolute_difference'):
                v.verify_solver_publication(path, *expected)

    def test_solver_archive_subspace_outside_original_threshold_fails(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            expected = self.solver_archive_fixture(path)
            clusters = v.read_csv(path / 'solver_subspaces.csv')
            clusters[0]['minimum_mass_correlation'] = 1 - 2e-6
            clusters[0]['maximum_principal_angle_degrees'] = float(np.degrees(np.arccos(1 - 2e-6)))
            v.write_csv(path / 'solver_subspaces.csv', clusters)
            with self.assertRaisesRegex(ValueError, 'subspace exceeds unchanged threshold'):
                v.verify_solver_publication(path, *expected)

    def test_solver_archive_residual_and_gram_thresholds_are_unchanged(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            expected = self.solver_archive_fixture(path)
            rows = v.read_csv(path / 'solver_comparison.csv')
            rows[0]['lapack_relative_residual'] = 2e-7
            v.write_csv(path / 'solver_comparison.csv', rows)
            with self.assertRaisesRegex(ValueError, 'residual exceeds unchanged threshold'):
                v.verify_solver_publication(path, *expected)
            expected = self.solver_archive_fixture(path)
            report = json.loads((path / 'report_data.json').read_text())
            report['solver_comparison']['synthetic']['m_orthogonality_max_error'] = 2e-5
            (path / 'report_data.json').write_text(json.dumps(report))
            with self.assertRaisesRegex(ValueError, 'M-orthogonality exceeds unchanged threshold'):
                v.verify_solver_publication(path, *expected)

    def test_solver_archive_indices_status_and_summary_consistency_stay_strict(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            expected = self.solver_archive_fixture(path)
            rows = v.read_csv(path / 'solver_comparison.csv')
            rows[0]['index'] = 99
            v.write_csv(path / 'solver_comparison.csv', rows)
            with self.assertRaises(ValueError):
                v.verify_solver_publication(path, *expected)
            expected = self.solver_archive_fixture(path)
            report = json.loads((path / 'report_data.json').read_text())
            report['solver_comparison']['synthetic']['maximum_eigenvalue_absolute_difference'] = 0
            (path / 'report_data.json').write_text(json.dumps(report))
            with self.assertRaisesRegex(ValueError, 'maximum_eigenvalue_absolute_difference'):
                v.verify_solver_publication(path, *expected)

    def test_saved_derived_tables_are_checked_against_reconstruction(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'derived.csv'
            v.write_csv(path, [dict(mode='TM11', value=1., status='failed_boundary')])
            v.verify_derived_table(path, [dict(mode='TM11', value=1., status='failed_boundary')])
            with self.assertRaises(ValueError):
                v.verify_derived_table(path, [dict(mode='TM11', value=2., status='failed_boundary')])
            with self.assertRaises(ValueError):
                v.verify_derived_table(path, [dict(mode='TM11', value=1., status='supported_candidate')])


if __name__ == '__main__':
    unittest.main()
