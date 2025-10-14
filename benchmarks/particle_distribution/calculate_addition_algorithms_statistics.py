import pandas as pd
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from contrib.python.scripts.aspect_data import read_statistics

# ------------------------------------Oscillating velocity------------------------------------ #

# Random, oscillating velocity
random_df = read_statistics('output_addition/output-random/statistics')
random_oscillating_data = [
    'Random, oscillating V',
    random_df['Cell Score Standard Deviation: '].mean(),
    random_df['Average particle distribution score: '].mean(),
    random_df['Maximal particle distribution score: '].mean()
]

# Histogram, oscillating velocity
histogram_df = read_statistics('output_addition/output-histogram/statistics')
histogram_oscillating_data = [
    'Histogram, oscillating V',
    histogram_df['Cell Score Standard Deviation: '].mean(),
    histogram_df['Average particle distribution score: '].mean(),
    histogram_df['Maximal particle distribution score: '].mean()
]

# Gaussian kernel function, oscillating velocity
gaussian_df = read_statistics('output_addition/output-gaussian/statistics')
gaussian_oscillating_data = [
    'Gaussian, oscillating V',
    gaussian_df['Cell Score Standard Deviation: '].mean(),
    gaussian_df['Average particle distribution score: '].mean(),
    gaussian_df['Maximal particle distribution score: '].mean(),
]

# Cutoff-w1 kernel function, oscillating velocity
cutoffw1_df = read_statistics('output_addition/output-cutoff-w1/statistics')
cutoffw1_oscillating_data = [
    'Cutoff-w1, oscillating V',
    cutoffw1_df['Cell Score Standard Deviation: '].mean(),
    cutoffw1_df['Average particle distribution score: '].mean(),
    cutoffw1_df['Maximal particle distribution score: '].mean()
]

# Cutoff-c1 kernel function, oscillating velocity
cutoffc1_df = read_statistics('output_addition/output-cutoff-c1/statistics')
cutoffc1_oscillating_data = [
    'Cutoff-c1, oscillating V',
    cutoffc1_df['Cell Score Standard Deviation: '].mean(),
    cutoffc1_df['Average particle distribution score: '].mean(),
    cutoffc1_df['Maximal particle distribution score: '].mean()
]

# Uniform kernel function, oscillating velocity
uniform_df = read_statistics('output_addition/output-uniform/statistics')
uniform_oscillating_data = [
    'Uniform, oscillating V',
    uniform_df['Cell Score Standard Deviation: '].mean(),
    uniform_df['Average particle distribution score: '].mean(),
    uniform_df['Maximal particle distribution score: '].mean()
]

# Triangular kernel function, oscillating velocity
triangular_df = read_statistics('output_addition/output-triangular/statistics')
triangular_oscillating_data = [
    'Triangular, oscillating V',
    triangular_df['Cell Score Standard Deviation: '].mean(),
    triangular_df['Average particle distribution score: '].mean(),
    triangular_df['Maximal particle distribution score: '].mean()
]

# ------------------------------------Constant velocity------------------------------------ #

# Random, constant velocity
random_constant_df = read_statistics('output_addition/output-random-constant-velocity/statistics')
random_constant_data = [
    'Random, constant V',
    random_constant_df['Cell Score Standard Deviation: '].mean(),
    random_constant_df['Average particle distribution score: '].mean(),
    random_constant_df['Maximal particle distribution score: '].mean()
]

# Histogram, constant velocity
histogram_constant_df = read_statistics('output_addition/output-histogram-constant-velocity/statistics')
histogram_constant_data = [
    'Histogram, constant V',
    histogram_constant_df['Cell Score Standard Deviation: '].mean(),
    histogram_constant_df['Average particle distribution score: '].mean(),
    histogram_constant_df['Maximal particle distribution score: '].mean()
]


# Gaussian kernel function, constant velocity
gaussian_constant_df = read_statistics('output_addition/output-gaussian-constant-velocity/statistics')
gaussian_constant_data = [
    'Gaussian, constant V',
    gaussian_constant_df['Cell Score Standard Deviation: '].mean(),
    gaussian_constant_df['Average particle distribution score: '].mean(),
    gaussian_constant_df['Maximal particle distribution score: '].mean()
]

# Cutoff-w1 kernel function, constant velocity
cutoffw1_constant_df = read_statistics('output_addition/output-cutoff-w1-constant-velocity/statistics')
cutoffw1_constant_data = [
    'Cutoff-w1, constant V',
    cutoffw1_constant_df['Cell Score Standard Deviation: '].mean(),
    cutoffw1_constant_df['Average particle distribution score: '].mean(),
    cutoffw1_constant_df['Maximal particle distribution score: '].mean(),
]

# Cutoff-c1 kernel function, constant velocity
cutoffc1_constant_df = read_statistics('output_addition/output-cutoff-c1-constant-velocity/statistics')
cutoffc1_constant_data = [
    'Cutoff-c1, constant V',
    cutoffc1_constant_df['Cell Score Standard Deviation: '].mean(),
    cutoffc1_constant_df['Average particle distribution score: '].mean(),
    cutoffc1_constant_df['Maximal particle distribution score: '].mean()
]

# Uniform kernel function, constant velocity
uniform_constant_df = read_statistics('output_addition/output-uniform-constant-velocity/statistics')
uniform_constant_data = [
    'Uniform, constant V',
    uniform_constant_df['Cell Score Standard Deviation: '].mean(),
    uniform_constant_df['Average particle distribution score: '].mean(),
    uniform_constant_df['Maximal particle distribution score: '].mean(),
]

# Triangular kernel function, constant velocity
triangular_constant_df = read_statistics('output_addition/output-triangular-constant-velocity/statistics')
triangular_constant_data = [
    'Triangular, constant V',
    triangular_constant_df['Cell Score Standard Deviation: '].mean(),
    triangular_constant_df['Average particle distribution score: '].mean(),
    triangular_constant_df['Maximal particle distribution score: '].mean()
]

output_data_array_oscillating = [
    random_oscillating_data,
    histogram_oscillating_data,
    gaussian_oscillating_data,
    cutoffw1_oscillating_data,
    cutoffc1_oscillating_data,
    uniform_oscillating_data,
    triangular_oscillating_data
]

output_data_array_constant = [
    random_constant_data,
    histogram_constant_data,
    gaussian_constant_data,
    cutoffw1_constant_data,
    cutoffc1_constant_data,
    uniform_constant_data,
    triangular_constant_data
]

column_names_output = [
    'Particle Addition Algorithm',
    'Time Averaged Cell Score Standard Deviation',
    'Time Averaged Mean Score',
    'Time Averaged Maximum Score',
]

output_dataframe_oscillating = pd.DataFrame(output_data_array_oscillating,columns=column_names_output)
output_dataframe_oscillating.to_csv('addition_algorithm_comparison_oscillatingV_data.csv',index=False)

output_dataframe_constant = pd.DataFrame(output_data_array_constant,columns=column_names_output)
output_dataframe_constant.to_csv('addition_algorithm_comparison_constantV_data.csv',index=False)
