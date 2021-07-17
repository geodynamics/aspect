import numpy as np

from aspect_postprocess_utils.errors import DataError

class AspectStateQuadrature:
    diff_tolerance = 1e-7

    def __init__(self, dimension, var_names, points, weights, data, time=0.0, diameter=None, q_type=None):
        self.dimension = dimension
        self.time = time
        self.diameter = diameter
        self.q_type = q_type
        self.var_names = var_names
        self.points = points
        self.weights = weights
        self.data = data

    def merge(self, other):
        """Merge datapoints from other file"""
        if not isinstance(other, AspectStateQuadrature):
            raise DataError("AspectStateQuadrature.merge requires AspectStateQuadrature argument")
        if not (self.dim == other.dim and
                self.time == other.time and
                self.diameter == other.diameter and
                self.q_type == other.q_type):
            raise DataError("Dataset metadata mismatch, files likely incorrect")
        for v1, v2 in zip(self.var_names, other.var_names):
            if v1 != v2:
                raise DataError(
                    "Dataset variable list mismatch {} and {}".format(
                        self.var_names,
                        v_names
                    )
                )
        self.points = np.concatenate((self.points, points), axis=0)
        self.weights = np.concatenate((self.points, points), axis=0)
        for v in self.var_names:
            self.data[v] = np.concatenate((self.data[v], data[v]), axis=0)

    def hash_order(self, hash_func):
        """Re-order dataset in order of ascending value of the provided function"""
        hashes = hash_func(self)

        self.reorder(hashes.argsort())

    def reorder(self, new_order):
        """Re-order dataset in terms of the given indicies"""
        self.points = self.points[new_order, :]
        self.weights = self.weights[new_order]
        for v in self.var_names:
            if v == 'V':
                self.data[v] = self.data[v][new_order, :]
            else:
                self.data[v] = self.data[v][new_order]


def box_point_hash(state):
    point_x = np.amin(state.points, axis=0)
    extent_x = np.amax(state.points, axis=0) - point_x
    count_x = np.zeros((state.dimension,), dtype=np.int64)
    for i in range(state.dimension):
        count_x[i] = np.unique(state.points[:, i]).shape[0]

    bin_size = extent_x/(count_x-1)

    point_x -= 0.5*bin_size

    hashes = np.zeros((state.points.shape[0],), dtype=np.int64)
    for i in range(state.dimension):
        hashes += (np.prod(count_x[:i])*np.floor((state.points[:,i]-point_x[i])/bin_size[i])).astype(np.int64)

    return hashes
