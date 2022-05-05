import numpy as np

from aspect_postprocess_utils.errors import DataError
from aspect_postprocess_utils.data import AspectStateQuadrature

class AspectLoad:
    @staticmethod
    def load_file(filename):
        raise NotImplementedError("No implementation provided for loading file in {}".format(__class__.__name__))

    @classmethod
    def load_files(cls, filenames):
        return_object = None
        for fn in filenames:
            if return_object is None:
                return_object = cls.load_file(fn)
            else:
                new_object = cls.load_file(fn)
                return_object.merge(new_object)
        return return_object


class AspectLoadASCII(AspectLoad):
    @staticmethod
    def load_file(filename):
        with open(filename, 'r') as f:
            k,v = f.readline().split(':')
            if(k == "DIMENSION"):
                dim = int(v)
            else:
                raise DataError("Invalid header of file {}".format(filename))

            k,v = f.readline().split(':')
            if(k == "TIME"):
                time = float(v)
            else:
                raise DataError("Invalid header of file {}".format(filename))

            k,v = f.readline().split(':')
            if(k == "CELL DIAMETER"):
                diameter = float(v)
            else:
                raise DataError("Invalid header of file {}".format(filename))

            k,v1,v2 = f.readline().split(':')
            if(k == "QUADRATURE"):
                q_type = (v1,int(v2))
            else:
                raise DataError("Invalid header of file {}".format(filename))

            head_vars = [v.strip() for v in f.readline().split('\t')]
            data_array = np.loadtxt(f, delimiter='\t')

            # Check dimension of point and velocity vectors
            if(not dim == len([v for v in head_vars if v.startswith('X_')])):
                raise DataError("Invalid header of file {}".format(filename))
            if(not dim == len([v for v in head_vars if v.startswith('V_')])):
                raise DataError("Invalid header of file {}".format(filename))

        head_vars = ['V'] + head_vars[dim*2+1:]
        data = {}
        points = data_array[:, 0:dim] # points
        weights = data_array[:, dim] # weights
        data['V'] = data_array[:, dim+1:2*dim+1] # velocities
        for i,v in enumerate(head_vars[1:]):
            data[v] = data_array[:,2*dim+1+i]
        return AspectStateQuadrature(dim, head_vars, points, weights, data,
                                     time=time, diameter=diameter, q_type=q_type)
