import csv
import numpy as np

from aspect_postprocess_utils.errors import DataError
from aspect_postprocess_utils.data import AspectStateQuadrature

class AspectStateDiff:
    def __init__(self, state_a, state_b):
        if not (isinstance(state_a, AspectStateQuadrature) and
                isinstance(state_b, AspectStateQuadrature)):
            raise DataError("Arguments to {this_class} constructor must be {data_class}".format(
                this_class=__class__.__name__,
                data_class=AspectStateQuadrature.__class__.__name__
            ))

        self.state_a = state_a
        self.state_b = state_b

        if self.state_a.weights.shape[0] != self.state_b.weights.shape[0]:
            raise DataError("Datasets not the same length".format(var))

        if self.diameter is None:
            diameter_tol = AspectStateQuadrature.diff_tolerance
        else:
            diameter_tol = AspectStateQuadrature.diff_tolerance * self.diameter

        point_diff_max = self.point_diff_ls()
        if not point_diff_max < diameter_tol:
            raise DataError("Dataset point order mismatch {diff} > {tol}".format(
                diff = point_diff_max,
                tol = diameter_tol,
            ))

    @property
    def variables(self):
        """List of variables included in both datasets"""
        return (set(self.state_a.var_names) &
                set(self.state_b.var_names))

    @property
    def diameter(self):
        """Minimum diameter of the two input state datasets"""
        if self.state_a.diameter is None or self.state_b.diameter is None:
            return None
        return np.minimum(self.state_a.diameter, self.state_b.diameter)

    def point_diff_ls(self):
        """Largest distance between the corresponding points in the two datasets."""
        return np.amax(np.linalg.norm(self.state_a.points-self.state_b.points, axis=1))

    def diff_l1(self, var):
        """Compute the L1 norm of the difference in the given variable between
        the datasets by quadrature."""
        if var == 'V':
            return np.sum(
                np.linalg.norm(self.state_a.data[var]-self.state_b.data[var], axis=1)
                *(self.state_a.weights+self.state_b.weights)/2.
            )
        else:
            return np.sum(
                np.abs(self.state_a.data[var]-self.state_b.data[var])
                *(self.state_a.weights+self.state_b.weights)/2.
            )

    def diff_l2(self, var):
        """Compute the L2 norm of the difference in the given variable between
        the datasets by quadrature."""
        if var == 'V':
            return np.sqrt(
                np.sum(
                    np.linalg.norm(self.state_a.data[var]-self.state_b.data[var], axis=1)**2
                    *(self.state_a.weights+self.state_b.weights)/2.
                )
            )
        else:
            return np.sqrt(
                np.sum(
                    np.abs(self.state_a.data[var]-self.state_b.data[var])**2
                    *(self.state_a.weights+self.state_b.weights)/2.
                )
            )

class StateConvergence:
    def __init__(self, variables=None, L1=True, L2=True):
        self.L1 = L1
        self.L2 = L2
        self.diff_vars = variables
        self.diff_data = []

    @property
    def diffs(self):
        """List of names for diff norms in dataset output"""
        return sorted(
            [n+"_L1" for n in self.diff_vars if self.L1] +
            [n+"_L2" for n in self.diff_vars if self.L2]
             )

    @property
    def rates(self):
        """List of names for diff norm convergence rates in dataset output"""
        return sorted(
            [n+"_L1_rate" for n in self.diff_vars if self.L1] +
            [n+"_L2_rate" for n in self.diff_vars if self.L2]
             )

    def add_diff(self, diff, label=""):
        # (AspectStateDiff, str) -> None
        """Add data from diff to convergence computation"""

        # If variables not specified, add all to list
        if self.diff_vars is None:
            self.diff_vars = diff.variables

        missing_variables = diff.variables-set(self.diff_vars)
        if missing_variables:
            raise DataError("Diff missing variables {}".format(missing_variables))

        row = {"label":label, "diameter": diff.diameter}
        for var in self.diff_vars:
            if self.L1:
                row[var+"_L1"] = diff.diff_l1(var)
            if self.L2:
                row[var+"_L2"] = diff.diff_l2(var)
        self.diff_data.append(row)

    def convergence_data(self, compute_rates=True):
        def diameter_key(d):
            return d['diameter']
        old_row = None
        for row in sorted(self.diff_data, key=diameter_key, reverse=True):
            out_row = row.copy()
            if old_row is None:
                for v in self.diffs:
                    out_row[v+"_rate"] = np.nan
            else:
                o_diam = old_row['diameter']
                n_diam = row['diameter']
                if o_diam is None or n_diam is None:
                    factor = 1.0
                else:
                    factor = np.log2(o_diam/n_diam)
                for v in self.diffs:
                    o_error = old_row[v]
                    n_error = row[v]
                    out_row[v+"_rate"] = (
                        np.log2(o_error/n_error)/factor
                    )
            old_row = row
            yield out_row

    def csv_convergence(self, ostream):
        fields = ['i', 'label'] + sorted(self.diffs + self.rates)
        writer = csv.DictWriter(ostream, fields)
        writer.writeheader()
        for i, row in enumerate(self.convergence_data()):
            write_row = {
                'i': i,
                'label': row['label'],
            }
            for k in self.diffs:
                write_row[k] = "{:.2e}".format(row[k])
            for k in self.rates:
                write_row[k] = "{:.2f}".format(row[k])
            writer.writerow(write_row)
