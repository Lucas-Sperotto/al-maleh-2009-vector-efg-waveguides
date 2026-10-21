#!/usr/bin/env python3
"""Reopen complete FEM/EFG exports, check algebra, and compare physical fields.

Requires numpy, scipy, matplotlib, and the system liblapacke shared library.
No execution of either discretization is hidden in this postprocessing step.
Run with --campaign DIR; --check-only reopens and verifies saved artifacts and
checks existing derived CSV/JSON without changing them.
"""
from __future__ import annotations

import argparse
import csv
import ctypes
import ctypes.util
import gzip
import itertools
import json
from pathlib import Path
import platform
import re
import sys

import numpy as np
from scipy import io, linalg, optimize, sparse
from scipy.sparse import csr_matrix

MODES = ((1, 1), (2, 1), (3, 1), (1, 2), (2, 2), (4, 1))
SINGLES = MODES[:4]
THRESHOLDS = dict(eigenvalue_relative_error=0.15, minimum_overlap=0.75,
                  divergence_ratio=0.25, normal_ratio=0.25,
                  algebraic_relative_residual=1e-7)


def read_csv(path):
    with Path(path).open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError('Refusing empty comparison table: ' + str(path))
    with Path(path).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_matrix(path):
    path = Path(path)
    with (gzip.open(path, 'rb') if path.suffix == '.gz' else path.open('rb')) as stream:
        matrix = io.mmread(stream)
    return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)


def finite(name, *arrays):
    for array in arrays:
        if not np.isfinite(np.asarray(array, dtype=float)).all():
            raise ValueError(name + ': nonfinite numerical artifact')


def eigenvalue(m, n):
    return np.pi ** 2 * (m * m / 4.0 + n * n)


def analytic_fields(points, modes=MODES):
    return np.column_stack([np.sin(m * np.pi * points[:, 0] / 2) *
                            np.sin(n * np.pi * points[:, 1]) for m, n in modes])


def normalize_fields(fields, weights):
    fields = np.asarray(fields, dtype=float)
    if fields.ndim == 1:
        fields = fields[:, None]
    norms = np.sqrt(np.sum(weights[:, None] * fields ** 2, axis=0))
    if np.any(norms <= np.finfo(float).tiny):
        raise ValueError('Cannot normalize a zero field')
    return fields / norms


def weighted_basis(fields, weights):
    """Orthonormal field basis in the common physical quadrature norm."""
    values = np.sqrt(weights)[:, None] * np.asarray(fields)
    u, singular, _ = linalg.svd(values, full_matrices=False)
    cutoff = max(values.shape) * np.finfo(float).eps * singular[0]
    if np.count_nonzero(singular > cutoff) != values.shape[1]:
        raise ValueError('Rank-deficient candidate field subspace')
    return u[:, :values.shape[1]]


def principal_metrics(left, right, weights):
    qleft, qright = weighted_basis(left, weights), weighted_basis(right, weights)
    singular = np.clip(linalg.svdvals(qleft.T @ qright), 0, 1)
    angles = np.degrees(np.arccos(singular))
    return dict(minimum_correlation=float(singular.min()),
                maximum_angle_degrees=float(angles.max()),
                projector_distance_frobenius=float(np.sqrt(max(
                    0, qleft.shape[1] + qright.shape[1] - 2 * np.sum(singular ** 2)))),
                singular_values=singular.tolist(), angles_degrees=angles.tolist())


def assign_modes(values, fields, points, weights, cutoff=0.0):
    """Assign simple modes by global overlap and the doublet by its span.

    A 50% spectral window limits exploratory candidates. It is wider than the
    unchanged 15% acceptance criterion. The eigenvalue never sets the assignment
    score for simple modes. The doublet maximizes its worst principal overlap.
    """
    references = normalize_fields(analytic_fields(points), weights)
    norms = np.sqrt(np.sum(weights[:, None] * fields ** 2, axis=0))
    valid = (values > cutoff) & (norms > 1e-12)
    normalized = np.zeros_like(fields)
    normalized[:, valid] = fields[:, valid] / norms[valid]
    overlaps = np.abs(references.T @ (weights[:, None] * normalized))
    costs = np.full((4, len(values)), 1e6)
    for row, mode in enumerate(SINGLES):
        allowed = valid & (np.abs(values / eigenvalue(*mode) - 1) <= 0.5)
        costs[row, allowed] = -overlaps[row, allowed]
    rows, columns = optimize.linear_sum_assignment(costs)
    if len(rows) != 4 or np.any(costs[rows, columns] >= 1e5):
        raise ValueError('No distinct field-based candidates for all simple modes')
    result = {MODES[row]: [int(column)] for row, column in zip(rows, columns)}
    target = eigenvalue(2, 2)
    candidates = [int(i) for i in np.flatnonzero(
        valid & (np.abs(values / target - 1) <= 0.5)) if i not in columns]
    if len(candidates) < 2:
        raise ValueError('No two-dimensional candidate span for TM22/TM41')
    qref = weighted_basis(references[:, 4:6], weights)
    best_pair, best_score = None, (-1.0, -np.inf)
    for pair in itertools.combinations(candidates, 2):
        try:
            qpair = weighted_basis(normalized[:, pair], weights)
        except ValueError:
            continue
        singular = linalg.svdvals(qref.T @ qpair)
        score = (float(singular.min()), -float(np.max(np.abs(values[list(pair)] / target - 1))))
        if score > best_score:
            best_pair, best_score = pair, score
    if best_pair is None:
        raise ValueError('No full-rank candidate doublet')
    result[(2, 2)] = result[(4, 1)] = list(best_pair)
    return result, normalized, overlaps


def p1_interpolation(nodes, elements, dof_nodes, points):
    """Locate points in the exported triangles, preserving mesh connectivity."""
    nodes, elements = np.asarray(nodes), np.asarray(elements, dtype=int)
    dof_nodes = np.asarray(dof_nodes, dtype=int)
    nodal_to_dof = np.full(len(nodes), -1, dtype=int)
    nodal_to_dof[dof_nodes] = np.arange(len(dof_nodes))
    assigned = np.zeros(len(points), dtype=bool)
    rows, columns, data = [], [], []
    for triangle in elements:
        vertices = nodes[triangle]
        low, high = vertices.min(axis=0) - 1e-12, vertices.max(axis=0) + 1e-12
        candidates = np.flatnonzero((~assigned) & np.all(points >= low, axis=1) &
                                   np.all(points <= high, axis=1))
        if not len(candidates):
            continue
        transform = np.column_stack((vertices[1] - vertices[0], vertices[2] - vertices[0]))
        uv = np.linalg.solve(transform, (points[candidates] - vertices[0]).T).T
        barycentric = np.column_stack((1 - uv[:, 0] - uv[:, 1], uv))
        inside = np.all(barycentric >= -1e-11, axis=1)
        for point, bary in zip(candidates[inside], barycentric[inside]):
            assigned[point] = True
            for node, coefficient in zip(triangle, bary):
                dof = nodal_to_dof[node]
                if dof >= 0:
                    rows.append(int(point)); columns.append(int(dof)); data.append(float(coefficient))
    if not assigned.all():
        raise ValueError('Common grid has points outside exported FEM mesh')
    return csr_matrix((data, (rows, columns)), shape=(len(points), len(dof_nodes)))


def pair_residuals(k, m, values, vectors):
    numerator = np.linalg.norm(k @ vectors - (m @ vectors) * values[None, :], axis=0)
    denominator = ((np.linalg.norm(k) + np.abs(values) * np.linalg.norm(m)) *
                   np.linalg.norm(vectors, axis=0))
    return numerator / np.maximum(denominator, 1e-300)


def audit_problem(path, family):
    values_table = read_csv(path / 'eigenvalues.csv')
    values = np.array([float(row['lambda']) for row in values_table])
    indices = [int(row['index']) for row in values_table]
    if indices != list(range(len(values))):
        raise ValueError(str(path) + ': spectrum index and vector-column convention broken')
    k = read_matrix(path / 'stiffness_matrix.mtx.gz')
    m = read_matrix(path / 'mass_matrix.mtx.gz')
    vectors = read_matrix(path / 'eigenvectors.mtx.gz')
    finite(str(path), values, k, m, vectors)
    if k.shape != m.shape or k.shape != vectors.shape or k.shape != (len(values), len(values)):
        raise ValueError(str(path) + ': complete spectrum/matrix dimensions do not agree')
    symmetry = {name: float(np.linalg.norm(mat - mat.T) / max(np.linalg.norm(mat), 1e-300))
                for name, mat in [('stiffness', k), ('mass', m)]}
    if max(symmetry.values()) > 1e-10:
        raise ValueError(str(path) + ': matrices are not symmetric')
    linalg.cholesky(m, lower=True)
    residuals = pair_residuals(k, m, values, vectors)
    gram_error = float(np.max(np.abs(vectors.T @ m @ vectors - np.eye(len(values)))))
    if residuals.max() > 1e-7 or gram_error > 1e-5:
        raise ValueError(str(path) + ': persisted eigenvectors fail residual or M-orthogonality')
    if np.any(np.diff(values) < -1e-10):
        raise ValueError(str(path) + ': unsorted complete eigenvalues')
    metadata = json.loads((path / 'run_metadata.json').read_text())
    geometry = metadata.get('geometry', {})
    if not isinstance(geometry, dict) or not {'a', 'b'} <= set(geometry):
        raise ValueError(str(path) + ': physical geometry a,b must be explicitly recorded')
    if family == 'fem' and (metadata.get('family') != 'TM' or
                            metadata.get('boundary_condition') != 'homogeneous Dirichlet Ez=0'):
        raise ValueError(str(path) + ': FEM artifact is not the TM Dirichlet family')
    if abs(float(geometry['a']) - 2) > 1e-12 or abs(float(geometry['b']) - 1) > 1e-12:
        raise ValueError(str(path) + ': physical geometry is not 2 x 1')
    cutoff = float(metadata.get('spectrum', {}).get('zero_cutoff', max(1e-9, values[-1] * 1e-9)))
    if family == 'efg':
        dofs = read_csv(path / 'dof_map.csv')
        if len(dofs) != len(values):
            raise ValueError('EFG degree-of-freedom map does not cover all coefficients')
        half = len(values) // 2
        if len(values) % 2 or [int(row['dof']) for row in dofs] != indices or \
                [row['component'] for row in dofs] != ['Hx'] * half + ['Hy'] * half or \
                [int(row['node']) for row in dofs] != list(range(half)) * 2:
            raise ValueError('EFG DOF convention must be all Hx then all Hy')
        positions = np.array([[float(row['x']), float(row['y']), float(row['support_radius'])] for row in dofs])
        finite('EFG degree-of-freedom map', positions)
        if not np.array_equal(positions[:half], positions[half:]) or (positions[:, 2] <= 0).any():
            raise ValueError('EFG component geometry/support maps differ')
        if not np.allclose(positions[:, :2].min(axis=0), [0, 0]) or not np.allclose(positions[:, :2].max(axis=0), [2, 1]):
            raise ValueError('EFG nodal geometry does not span the 2 x 1 rectangle')
        diagnostics = read_csv(path / 'mode_diagnostics.csv')
        if [int(row['index']) for row in diagnostics] != indices:
            raise ValueError('EFG diagnostics do not follow eigenvector columns')
        for row in diagnostics:
            finite('EFG persisted diagnostics', [float(value) for key, value in row.items()
                   if key != 'status' and value != ''])
            if int(row['index']) >= len(values) or float(row['lambda']) != values[int(row['index'])]:
                raise ValueError('EFG diagnostic eigenvalue differs from its vector column')
    else:
        diagnostics = []
    summary = dict(path=str(path), dofs=len(values), complete_vectors=vectors.shape[1],
                   maximum_relative_residual=float(residuals.max()),
                   m_orthogonality_max_error=gram_error, symmetry_relative_frobenius=symmetry,
                   minimum_mass_eigenvalue=float(linalg.eigvalsh(m, subset_by_index=[0, 0])[0]),
                   zero_cutoff=cutoff, numerical_kernel_dimension=int(np.sum(values <= cutoff)),
                   retained_nonpositive_eigenvalues=int(np.sum(values <= 0)))
    return dict(path=path, values=values, vectors=vectors, k=k, m=m,
                residuals=residuals, metadata=metadata, cutoff=cutoff,
                diagnostics=diagnostics, audit=summary)


def lapack_dsygv(k, m):
    """Independent LAPACKE dsygv, double precision, exact persisted matrices."""
    library = ctypes.util.find_library('lapacke')
    if not library:
        raise RuntimeError('liblapacke is required for the independent EFG solver check')
    api = ctypes.CDLL(library).LAPACKE_dsygv
    pointer = ctypes.POINTER(ctypes.c_double)
    api.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_char, ctypes.c_char,
                    ctypes.c_int, pointer, ctypes.c_int, pointer, ctypes.c_int, pointer]
    api.restype = ctypes.c_int
    a, b = np.array(k, order='F', copy=True), np.array(m, order='F', copy=True)
    values = np.empty(len(k), dtype=np.float64)
    status = api(102, 1, b'V', b'U', len(k), a.ctypes.data_as(pointer), len(k),
                 b.ctypes.data_as(pointer), len(k), values.ctypes.data_as(pointer))
    if status != 0:
        raise RuntimeError('LAPACKE_dsygv failed with info=' + str(status))
    return values, a, library


def solver_comparison(problem, label):
    # Eigen's reference solver performs precisely these averages before solve.
    stiffness = .5 * (problem['k'] + problem['k'].T)
    mass = .5 * (problem['m'] + problem['m'].T)
    values, vectors, library = lapack_dsygv(stiffness, mass)
    original = problem['values']
    residuals = pair_residuals(stiffness, mass, values, vectors)
    gram_error = float(np.max(np.abs(vectors.T @ mass @ vectors - np.eye(len(values)))))
    # Backward-error scale follows the whitened self-adjoint operator norm.
    # Absolute comparison is mandatory near zero; relative errors there are meaningless.
    spectrum_scale = max(float(np.max(np.abs(original))), float(np.max(np.abs(values))), 1.0)
    absolute_tolerance = 200 * len(values) * np.finfo(float).eps * spectrum_scale
    tolerance = absolute_tolerance + 1e-9 * np.maximum(np.abs(original), np.abs(values))
    if np.any(np.abs(values - original) > tolerance):
        raise ValueError(label + ': independent LAPACK eigenvalues exceed absolute+relative tolerance')
    if residuals.max() > 1e-7 or gram_error > 1e-5:
        raise ValueError(label + ': independent LAPACK residual/orthogonality failure')
    rows = [dict(case=label, index=i, eigen_lambda=float(original[i]),
                 lapack_lambda=float(values[i]), absolute_difference=float(abs(original[i] - values[i])),
                 relative_difference=float(abs(original[i] - values[i]) / max(abs(original[i]), abs(values[i]), problem['cutoff'])),
                 allowed_absolute_difference=float(tolerance[i]),
                 classification='numerical_kernel' if original[i] <= problem['cutoff'] else 'positive',
                 eigen_relative_residual=float(problem['residuals'][i]),
                 lapack_relative_residual=float(residuals[i]), status='pass') for i in range(len(values))]
    # Treat the entire numerical kernel as a subspace; individual basis vectors
    # in it and in repeated positive clusters are arbitrary.
    groups = []
    kernel = np.flatnonzero(original <= problem['cutoff'])
    if len(kernel):
        groups.append(('numerical_kernel', kernel))
    positive = np.flatnonzero(original > problem['cutoff'])
    group = []
    for i in positive:
        if group and abs(original[i] - original[group[-1]]) > absolute_tolerance + 1e-8 * max(1, abs(original[i])):
            groups.append(('positive_cluster', np.array(group))); group = []
        group.append(i)
    if group:
        groups.append(('positive_cluster', np.array(group)))
    clusters = []
    for kind, group in groups:
        singular = np.clip(linalg.svdvals(problem['vectors'][:, group].T @ mass @ vectors[:, group]), 0, 1)
        min_correlation = float(singular.min())
        if min_correlation < 1 - 1e-6:
            raise ValueError(label + ': LAPACK/Eigen invariant-subspace mismatch')
        clusters.append(dict(case=label, classification=kind, first_index=int(group[0]),
                             last_index=int(group[-1]), dimension=len(group),
                             minimum_mass_correlation=min_correlation,
                             maximum_principal_angle_degrees=float(np.degrees(np.arccos(min_correlation))),
                             status='pass'))
    summary = dict(backend='LAPACKE_dsygv', library=library,
                   input='A=0.5*(persisted_K+persisted_K.T), B=0.5*(persisted_M+persisted_M.T); exactly the Eigen input transformation',
                   eigenvalue_absolute_tolerance=absolute_tolerance,
                   eigenvalue_relative_tolerance=1e-9,
                   maximum_eigenvalue_absolute_difference=float(np.max(np.abs(values - original))),
                   maximum_relative_residual=float(residuals.max()),
                   m_orthogonality_max_error=gram_error,
                   minimum_subspace_mass_correlation=min(r['minimum_mass_correlation'] for r in clusters))
    return rows, clusters, summary


def load_fields(problem, family, points, weights):
    path = problem['path']
    if family == 'efg':
        grid = read_csv(path / 'grid.csv')
        local_points = np.array([[float(row['x']), float(row['y'])] for row in grid])
        local_weights = np.array([float(row['weight']) for row in grid])
        if local_points.shape != points.shape or not np.allclose(local_points, points, atol=1e-14, rtol=0) or not np.allclose(local_weights, weights):
            raise ValueError('EFG cases do not share the same physical grid')
        fields = read_matrix(path / 'field_curl.mtx.gz')
    else:
        node_rows = read_csv(path / 'mesh_nodes.csv')
        nodes = np.array([[float(row['x']), float(row['y'])] for row in node_rows])
        element_rows = read_csv(path / 'mesh_elements.csv')
        elements = np.array([[int(row['node0']), int(row['node1']), int(row['node2'])] for row in element_rows])
        dof_rows = read_csv(path / 'dof_map.csv')
        dof_nodes = np.array([int(row['node_index']) for row in dof_rows])
        if len(dof_nodes) != len(problem['values']):
            raise ValueError('FEM free-DOF map is incomplete')
        if not np.allclose(nodes.min(axis=0), [0, 0]) or not np.allclose(nodes.max(axis=0), [2, 1]):
            raise ValueError('FEM exported mesh geometry differs from the EFG rectangle')
        areas = np.abs(np.cross(nodes[elements[:, 1]] - nodes[elements[:, 0]],
                                nodes[elements[:, 2]] - nodes[elements[:, 0]])) / 2
        if np.any(areas <= 0) or abs(areas.sum() - 2) > 1e-12:
            raise ValueError('FEM triangle areas do not cover the physical rectangle')
        interpolation = p1_interpolation(nodes, elements, dof_nodes, points)
        fields = interpolation @ problem['vectors']
        problem['audit']['mesh_nodes'] = len(nodes)
        problem['audit']['mesh_triangles'] = len(elements)
        problem['audit']['dirichlet_boundary_exact_zero'] = bool(np.all(
            (nodes[dof_nodes, 0] > 0) & (nodes[dof_nodes, 0] < 2) &
            (nodes[dof_nodes, 1] > 0) & (nodes[dof_nodes, 1] < 1)))
    finite('reconstructed ' + family, fields)
    if fields.shape != (len(points), len(problem['values'])):
        raise ValueError('Field columns do not follow full spectrum indices')
    assignments, normalized, overlaps = assign_modes(problem['values'], fields, points, weights, problem['cutoff'])
    problem.update(fields=fields, normalized=normalized, assignments=assignments, overlaps=overlaps)
    return problem


def mode_record(problem, mode, points, weights):
    indices = problem['assignments'][mode]
    values = problem['values'][indices]
    reference = analytic_fields(points, [mode])
    if len(indices) == 1:
        correlation = float(problem['overlaps'][MODES.index(mode), indices[0]])
        maximum_angle = float(np.degrees(np.arccos(np.clip(correlation, 0, 1))))
    else:
        metrics = principal_metrics(analytic_fields(points, MODES[4:]), problem['fields'][:, indices], weights)
        correlation = metrics['minimum_correlation']; maximum_angle = metrics['maximum_angle_degrees']
    record = dict(indices=';'.join(map(str, indices)), lambda_value=float(values.mean()),
                  lambda_min=float(values.min()), lambda_max=float(values.max()),
                  relative_error=float(np.max(np.abs(values / eigenvalue(*mode) - 1))),
                  analytic_overlap=correlation, analytic_maximum_angle_degrees=maximum_angle,
                  residual=float(problem['residuals'][indices].max()))
    if problem['diagnostics']:
        record['divergence_ratio'] = max(float(problem['diagnostics'][i]['divergence_ratio']) for i in indices)
        record['normal_ratio'] = max(float(problem['diagnostics'][i]['normal_ratio']) for i in indices)
        criteria = dict(eigenvalue=record['relative_error'] <= THRESHOLDS['eigenvalue_relative_error'],
                        field=correlation >= THRESHOLDS['minimum_overlap'],
                        divergence=record['divergence_ratio'] <= THRESHOLDS['divergence_ratio'],
                        boundary=record['normal_ratio'] <= THRESHOLDS['normal_ratio'],
                        residual=record['residual'] <= THRESHOLDS['algebraic_relative_residual'])
    else:
        record['divergence_ratio'] = ''; record['normal_ratio'] = ''
        criteria = dict(eigenvalue=record['relative_error'] <= THRESHOLDS['eigenvalue_relative_error'],
                        field=correlation >= THRESHOLDS['minimum_overlap'],
                        boundary=problem['audit']['dirichlet_boundary_exact_zero'],
                        residual=record['residual'] <= THRESHOLDS['algebraic_relative_residual'])
    failed = [name for name, passes in criteria.items() if not passes]
    record['status'] = 'supported_candidate' if not failed else 'failed_' + '_'.join(failed)
    record['matching'] = 'field_global_assignment' if len(indices) == 1 else 'joint_2d_subspace_centroid'
    return record


def compare(efg, fem, points, weights):
    rows = []
    for mode in MODES:
        er = mode_record(efg, mode, points, weights)
        fr = mode_record(fem, mode, points, weights)
        ei, fi = efg['assignments'][mode], fem['assignments'][mode]
        if len(ei) == 1:
            correlation = float(abs(np.sum(weights * efg['normalized'][:, ei[0]] * fem['normalized'][:, fi[0]])))
            metrics = dict(minimum_correlation=correlation,
                           maximum_angle_degrees=float(np.degrees(np.arccos(np.clip(correlation, 0, 1)))),
                           projector_distance_frobenius=float(np.sqrt(max(0, 2 - 2 * correlation ** 2))))
            field_error = float(np.sqrt(max(0, 2 - 2 * correlation)))
        else:
            metrics = principal_metrics(efg['fields'][:, ei], fem['fields'][:, fi], weights)
            field_error = float(np.sqrt(max(0, 2 - 2 * metrics['minimum_correlation'])))
        row = dict(mode='TM%d%d' % mode, m=mode[0], n=mode[1], analytic_lambda=eigenvalue(*mode))
        row.update({'efg_' + key: value for key, value in er.items()})
        row.update({'fem_' + key: value for key, value in fr.items()})
        row.update(efg_fem_relative_difference=abs(er['lambda_value'] - fr['lambda_value']) / eigenvalue(*mode),
                   field_correlation=metrics['minimum_correlation'],
                   sign_aligned_relative_l2_field_error=field_error,
                   field_maximum_principal_angle_degrees=metrics['maximum_angle_degrees'],
                   field_projector_distance_frobenius=metrics['projector_distance_frobenius'],
                   status='supported_cross_method_candidate' if
                   er['status'] == fr['status'] == 'supported_candidate' and
                   metrics['minimum_correlation'] >= THRESHOLDS['minimum_overlap']
                   else 'physical_validation_failed')
        rows.append(row)
    return rows


def make_figures(destination, convergence, comparisons, efg_runs, fem_runs, points, weights, sensitivity):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    destination.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 9, 'figure.dpi': 130, 'savefig.bbox': 'tight'})
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    for axis, mode in zip(axes.flat, MODES):
        name = 'TM%d%d' % mode
        for family in ('EFG', 'FEM'):
            rows = [r for r in convergence if r['mode'] == name and r['method'] == family]
            axis.plot([r['dofs'] for r in rows], [r['lambda_value'] for r in rows], 'o-', label=family)
        axis.axhline(eigenvalue(*mode), color='black', ls=':', label='Analytic')
        axis.set(title=name, xlabel='DOFs', ylabel=r'$\lambda$'); axis.legend()
    fig.suptitle('Independent eigenvalue convergence; doublet values are centroids')
    fig.tight_layout(); fig.savefig(destination / 'eigenvalue_convergence.png'); plt.close(fig)
    fig, axis = plt.subplots(figsize=(9, 4))
    x = np.arange(len(MODES))
    axis.bar(x - .18, [100 * row['efg_relative_error'] for row in comparisons], .36, label='EFG')
    axis.bar(x + .18, [100 * row['fem_relative_error'] for row in comparisons], .36, label='FEM')
    axis.axhline(15, color='black', ls=':', label='15% criterion')
    axis.set(xticks=x, xticklabels=[r['mode'] for r in comparisons], ylabel='Maximum eigenvalue error (%)')
    axis.legend(); fig.tight_layout(); fig.savefig(destination / 'eigenvalue_errors.png'); plt.close(fig)
    efg, fem = efg_runs[-1], fem_runs[-1]
    unique_x, unique_y = np.unique(points[:, 0]), np.unique(points[:, 1])
    fig, axes = plt.subplots(3, 3, figsize=(12, 7), sharex=True, sharey=True)
    for row, mode in enumerate(((1, 1), (2, 1), (1, 2))):
        reference = normalize_fields(analytic_fields(points, [mode]), weights)[:, 0]
        fields = [reference]
        for problem in (efg, fem):
            field = problem['normalized'][:, problem['assignments'][mode][0]].copy()
            field *= 1 if np.dot(weights * reference, field) >= 0 else -1
            fields.append(field)
        maximum = max(abs(field).max() for field in fields)
        for col, (title, field) in enumerate(zip(('Analytic Ez', 'EFG curl(H)', 'FEM P1 Ez'), fields)):
            # Export grid is y-major, x-fast. Sorting supports other equivalent row orders.
            order = np.lexsort((points[:, 0], points[:, 1]))
            image = axes[row, col].pcolormesh(unique_x, unique_y,
                field[order].reshape(len(unique_y), len(unique_x)), shading='nearest',
                cmap='RdBu_r', vmin=-maximum, vmax=maximum)
            axes[row, col].set(title=('TM%d%d — ' % mode) + title, aspect='equal', xlabel='x', ylabel='y')
            fig.colorbar(image, ax=axes[row, col], fraction=.035)
    fig.suptitle('L2-normalized fields on the shared grid; signs aligned to analytic Ez')
    fig.tight_layout(); fig.savefig(destination / 'field_comparison.png'); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for problem in efg_runs:
        records = [mode_record(problem, mode, points, weights) for mode in MODES]
        for axis, key in zip(axes, ('divergence_ratio', 'normal_ratio')):
            axis.plot(np.arange(6), [row[key] for row in records], 'o-', label=problem['path'].name)
            axis.axhline(.25, color='black', ls=':')
            axis.set(xticks=np.arange(6), xticklabels=['TM%d%d' % mode for mode in MODES], ylabel=key)
            axis.legend()
    fig.tight_layout(); fig.savefig(destination / 'efg_physical_indicators.png'); plt.close(fig)
    fig, axis = plt.subplots(figsize=(8, 4))
    for problem in efg_runs:
        rows = [r for r in sensitivity if r['case'] == problem['path'].name]
        axis.semilogx([r['zero_relative'] for r in rows], [r['kernel_dimension'] for r in rows], 'o-', label=problem['path'].name)
    axis.set(xlabel='Relative numerical-zero cutoff (max with absolute 1e-9)', ylabel='Classified kernel dimension')
    axis.legend(); fig.tight_layout(); fig.savefig(destination / 'kernel_threshold_sensitivity.png'); plt.close(fig)


def check_published_tables(campaign):
    count = 0
    for filename in ('comparison.csv', 'convergence.csv', 'solver_comparison.csv',
                     'solver_subspaces.csv', 'kernel_sensitivity.csv'):
        rows = read_csv(campaign / filename)
        if not rows:
            raise ValueError('Empty published table: ' + filename)
        for row in rows:
            for value in row.values():
                if value.strip().lower() in ('nan', 'inf', '-inf', '+inf'):
                    raise ValueError('Nonfinite published value in ' + filename)
        count += len(rows)
    report = json.loads((campaign / 'report_data.json').read_text(), parse_constant=lambda token:
                        (_ for _ in ()).throw(ValueError('Nonfinite JSON: ' + token)))
    if report['thresholds'] != THRESHOLDS:
        raise ValueError('Published acceptance criteria changed')
    return count


def verify_derived_table(path, expected, *, ignored_columns=()):
    """Detect stale or manually edited derived data, beyond a finiteness check."""
    actual = read_csv(path)
    if len(actual) != len(expected) or set(actual[0]) != set(expected[0]):
        raise ValueError('Published table schema/row count differs: ' + str(path))
    for row, reference in zip(actual, expected):
        for key, value in reference.items():
            if key in ignored_columns:
                continue
            if isinstance(value, (int, float, np.integer, np.floating)):
                absolute_tolerance = 2e-6 if 'angle_degrees' in key else 1e-12
                if not np.isclose(float(row[key]), value, rtol=1e-9, atol=absolute_tolerance):
                    raise ValueError('Published table differs from recomputation: ' + str(path) + ':' + key)
            elif row[key] != str(value):
                raise ValueError('Published classification differs from recomputation: ' + str(path) + ':' + key)


def verify_solver_publication(campaign, expected_rows, expected_clusters, expected_summaries, cutoffs):
    """Validate archived and newly solved LAPACK results against the same gates.

    A different BLAS implementation changes roundoff errors and near-zero
    principal angles. Archived measurements remain intact: they must satisfy
    the original per-eigenvalue tolerances and agree internally with their
    published summaries. Indices, Eigen data, classifications and pass/fail
    statuses retain strict comparison with reconstruction.
    """
    campaign = Path(campaign)
    rows_path = campaign / 'solver_comparison.csv'
    clusters_path = campaign / 'solver_subspaces.csv'
    verify_derived_table(rows_path, expected_rows, ignored_columns={
        'lapack_lambda', 'absolute_difference', 'relative_difference',
        'allowed_absolute_difference', 'lapack_relative_residual'})
    verify_derived_table(clusters_path, expected_clusters, ignored_columns={
        'minimum_mass_correlation', 'maximum_principal_angle_degrees'})
    rows, clusters = read_csv(rows_path), read_csv(clusters_path)
    archived = json.loads((campaign / 'report_data.json').read_text())['solver_comparison']
    if set(archived) != set(expected_summaries) or set(cutoffs) != set(expected_summaries):
        raise ValueError('Published solver-summary cases differ')

    def close(actual, expected, label):
        # Round-trip consistency of measurements from one archived execution,
        # not equality of errors produced by different LAPACK/BLAS builds.
        if not np.isfinite(actual) or not np.isclose(actual, expected, rtol=1e-9, atol=1e-15):
            raise ValueError('Published solver measurement is inconsistent: ' + label)

    for case, expected in expected_summaries.items():
        summary = archived[case]
        if set(summary) != set(expected):
            raise ValueError('Published solver-summary schema differs: ' + case)
        for key in ('backend', 'input'):
            if summary[key] != expected[key]:
                raise ValueError('Published solver backend/input differs: ' + case)
        if not isinstance(summary['library'], str) or not summary['library']:
            raise ValueError('Published solver library provenance is missing: ' + case)
        for key in ('eigenvalue_absolute_tolerance', 'eigenvalue_relative_tolerance'):
            close(float(summary[key]), expected[key], case + '.' + key)
        case_rows = [row for row in rows if row['case'] == case]
        case_clusters = [row for row in clusters if row['case'] == case]
        for row in case_rows:
            original, lapack = float(row['eigen_lambda']), float(row['lapack_lambda'])
            absolute_error = abs(original - lapack)
            allowed = (expected['eigenvalue_absolute_tolerance'] +
                       expected['eigenvalue_relative_tolerance'] * max(abs(original), abs(lapack)))
            if not np.isfinite(lapack) or absolute_error > allowed:
                raise ValueError('Published LAPACK eigenvalue exceeds declared tolerance: ' + case + ':' + row['index'])
            close(float(row['allowed_absolute_difference']), allowed, case + '.allowed_absolute_difference')
            close(float(row['absolute_difference']), absolute_error, case + '.absolute_difference')
            relative_error = absolute_error / max(abs(original), abs(lapack), cutoffs[case])
            close(float(row['relative_difference']), relative_error, case + '.relative_difference')
            residual = float(row['lapack_relative_residual'])
            if not np.isfinite(residual) or not 0 <= residual <= THRESHOLDS['algebraic_relative_residual']:
                raise ValueError('Published LAPACK residual exceeds unchanged threshold: ' + case)
        for row in case_clusters:
            correlation = float(row['minimum_mass_correlation'])
            if not np.isfinite(correlation) or not 1 - 1e-6 <= correlation <= 1:
                raise ValueError('Published LAPACK subspace exceeds unchanged threshold: ' + case)
            angle = float(np.degrees(np.arccos(correlation)))
            close(float(row['maximum_principal_angle_degrees']), angle, case + '.maximum_principal_angle_degrees')
        close(float(summary['maximum_eigenvalue_absolute_difference']),
              max(float(row['absolute_difference']) for row in case_rows), case + '.maximum_eigenvalue_absolute_difference')
        close(float(summary['maximum_relative_residual']),
              max(float(row['lapack_relative_residual']) for row in case_rows), case + '.maximum_relative_residual')
        close(float(summary['minimum_subspace_mass_correlation']),
              min(float(row['minimum_mass_correlation']) for row in case_clusters), case + '.minimum_subspace_mass_correlation')
        gram_error = float(summary['m_orthogonality_max_error'])
        if not np.isfinite(gram_error) or not 0 <= gram_error <= 1e-5:
            raise ValueError('Published LAPACK M-orthogonality exceeds unchanged threshold: ' + case)
    return archived


def verify_report(path, expected, *, solver_archive=None):
    """Verify published substantive JSON data; execution environment may differ."""
    actual = json.loads(Path(path).read_text())
    actual.pop('environment', None)
    expected = {key: value for key, value in expected.items() if key != 'environment'}
    if solver_archive is not None:
        # This section was independently checked against the original gates
        # and the archived per-pair/per-cluster tables, above.
        expected['solver_comparison'] = solver_archive

    def compare(actual_value, expected_value, location):
        if isinstance(expected_value, dict):
            if not isinstance(actual_value, dict) or set(actual_value) != set(expected_value):
                raise ValueError('Published report keys differ: ' + location)
            for key in expected_value:
                compare(actual_value[key], expected_value[key], location + '.' + key)
        elif isinstance(expected_value, list):
            if not isinstance(actual_value, list) or len(actual_value) != len(expected_value):
                raise ValueError('Published report lengths differ: ' + location)
            for index, (actual_item, expected_item) in enumerate(zip(actual_value, expected_value)):
                compare(actual_item, expected_item, location + '[' + str(index) + ']')
        elif isinstance(expected_value, (float, int)) and not isinstance(expected_value, bool):
            tolerance = 2e-6 if 'angle_degrees' in location else 1e-12
            if not isinstance(actual_value, (float, int)) or not np.isclose(actual_value, expected_value, rtol=1e-9, atol=tolerance):
                raise ValueError('Published report numerical value differs: ' + location)
        elif actual_value != expected_value:
            raise ValueError('Published report value differs: ' + location)
    compare(actual, expected, str(path))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--no-figures', action='store_true', help='Write tables and audits without rendering figures')
    args = parser.parse_args(argv)
    campaign = args.campaign.resolve()
    efg_paths = sorted((p for p in (campaign / 'efg').iterdir()
                        if p.is_dir() and re.fullmatch(r'\d+x\d+', p.name)),
                       key=lambda p: int(p.name.split('x')[0]))
    fem_paths = [campaign / 'fem' / name for name in ('coarse', 'medium', 'fine')]
    if len(efg_paths) != 3 or not all(path.is_dir() for path in fem_paths):
        raise ValueError('Expected three complete EFG and FEM refinements')
    grid = read_csv(efg_paths[0] / 'grid.csv')
    points = np.array([[float(row['x']), float(row['y'])] for row in grid])
    weights = np.array([float(row['weight']) for row in grid])
    finite('common physical grid', points, weights)
    if (weights <= 0).any() or abs(weights.sum() - 2) > 1e-12 or not np.all((points > [0, 0]) & (points < [2, 1])):
        raise ValueError('Common quadrature must cover interior of the 2 x 1 rectangle')
    efg_runs, fem_runs, solver_rows, cluster_rows, solver_summaries = [], [], [], [], {}
    for family, paths, destination in [('efg', efg_paths, efg_runs), ('fem', fem_paths, fem_runs)]:
        for path in paths:
            print('Reopening', family, path.name, flush=True)
            problem = load_fields(audit_problem(path, family), family, points, weights)
            destination.append(problem)
            if family == 'efg':
                rows, clusters, summary = solver_comparison(problem, path.name)
                solver_rows.extend(rows); cluster_rows.extend(clusters); solver_summaries[path.name] = summary
    convergence = []
    for family, runs in [('EFG', efg_runs), ('FEM', fem_runs)]:
        for problem in runs:
            for mode in MODES:
                record = mode_record(problem, mode, points, weights)
                convergence.append(dict(method=family, case=problem['path'].name,
                                        mode='TM%d%d' % mode, dofs=len(problem['values']),
                                        analytic_lambda=eigenvalue(*mode), **record))
    comparisons = compare(efg_runs[-1], fem_runs[-1], points, weights)
    sensitivity = []
    for problem in efg_runs:
        for relative in (1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4):
            cutoff = max(1e-9, relative * max(0, float(problem['values'].max())))
            sensitivity.append(dict(case=problem['path'].name, zero_absolute=1e-9,
                                    zero_relative=relative, effective_cutoff=cutoff,
                                    kernel_dimension=int(np.sum(problem['values'] <= cutoff)),
                                    positive_dimension=int(np.sum(problem['values'] > cutoff)),
                                    materially_negative_dimension=int(np.sum(problem['values'] < -cutoff)),
                                    changes_solution='false'))
    report = dict(thresholds=THRESHOLDS,
                  environment=dict(python=sys.version, numpy=np.__version__, platform=platform.platform()),
                  common_grid=dict(points=len(points), area=float(weights.sum()), field_normalization='weighted physical L2',
                                   fem_field='P1 longitudinal Ez, exact barycentric interpolation on exported triangles',
                                   efg_field='MLS transverse curl(H), evaluated by EFG executable'),
                  assignment=dict(simple='global maximum absolute analytic field overlap in 50% spectral window',
                                  doublet='joint two-dimensional span maximizing minimum analytic principal overlap',
                                  acceptance='unchanged M4 thresholds; doublet eigenvalue reports centroid and min/max'),
                  audits={family: [dict(problem['audit'], path=str(problem['path'].relative_to(campaign))) for problem in runs]
                          for family, runs in [('efg', efg_runs), ('fem', fem_runs)]},
                  solver_comparison=solver_summaries,
                  residual_definition='||Kv-lambda Mv||2 / ((||K||F+abs(lambda)||M||F)||v||2)',
                  comparisons=comparisons)
    if args.check_only:
        count = check_published_tables(campaign)
        solver_archive = verify_solver_publication(
            campaign, solver_rows, cluster_rows, solver_summaries,
            {problem['path'].name: problem['cutoff'] for problem in efg_runs})
        verify_report(campaign / 'report_data.json', report, solver_archive=solver_archive)
        verify_derived_table(campaign / 'comparison.csv', comparisons)
        verify_derived_table(campaign / 'convergence.csv', convergence)
        verify_derived_table(campaign / 'kernel_sensitivity.csv', sensitivity)
        for problem in fem_runs:
            stored = read_matrix(problem['path'] / 'common_grid_ez.mtx.gz')
            columns = read_csv(problem['path'] / 'common_grid_columns.csv')
            indices = [int(row['eigenvalue_index']) for row in columns]
            if not np.allclose(stored, problem['fields'][:, indices], rtol=1e-13, atol=1e-13):
                raise ValueError('Persisted FEM common-grid fields do not reconstruct from complete eigenvectors')
        print('PASS: 6 complete exports, 3 independent LAPACKE solves, field reconstruction, assignment, and', count, 'published rows')
        return 0
    write_csv(campaign / 'comparison.csv', comparisons)
    write_csv(campaign / 'convergence.csv', convergence)
    write_csv(campaign / 'solver_comparison.csv', solver_rows)
    write_csv(campaign / 'solver_subspaces.csv', cluster_rows)
    write_csv(campaign / 'kernel_sensitivity.csv', sensitivity)
    (campaign / 'report_data.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    for problem in fem_runs:
        indices = sorted(set(i for mode_indices in problem['assignments'].values() for i in mode_indices))
        with gzip.open(problem['path'] / 'common_grid_ez.mtx.gz', 'wb') as stream:
            io.mmwrite(stream, problem['fields'][:, indices], precision=17)
        write_csv(problem['path'] / 'common_grid_columns.csv',
                  [dict(field_column=i, eigenvalue_index=index) for i, index in enumerate(indices)])
    if not args.no_figures:
        make_figures(campaign / 'figures', convergence, comparisons, efg_runs, fem_runs, points, weights, sensitivity)
    print('Wrote comparisons, convergence, solver checks, sensitivity, report_data.json and 5 figures to', campaign)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
