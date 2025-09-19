import pandas as pd
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from contrib.python.scripts.aspect_data import read_statistics

# Random selection / no kernel function, oscillating velocity
#random_df = pd.read_table('output-random/statistics',names=column_names,comment='#',engine='python',sep='\s+')
random_df = read_statistics('output-random/statistics')
random_oscillating_data = [
    'Random, oscillating V',
    random_df['Cell Score Standard Deviation: '].mean(),
    random_df['Average particle distribution score: '].mean(),
    random_df['Maximal particle distribution score: '].mean()
]

# Gaussian kernel function, oscillating velocity
gaussian_df = read_statistics('output-gaussian/statistics')
gaussian_oscillating_data = [
    'Gaussian, oscillating V',
    gaussian_df['Cell Score Standard Deviation: '].mean(),
    gaussian_df['Average particle distribution score: '].mean(),
    gaussian_df['Maximal particle distribution score: '].mean(),
]

# Cutoff-w1 kernel function, oscillating velocity
cutoffw1_df = read_statistics('output-cutoff-w1/statistics')
cutoffw1_oscillating_data = [
    'Cutoff-w1, oscillating V',
    cutoffw1_df['Cell Score Standard Deviation: '].mean(),
    cutoffw1_df['Average particle distribution score: '].mean(),
    cutoffw1_df['Maximal particle distribution score: '].mean()
]

# Cutoff-c1 kernel function, oscillating velocity
cutoffc1_df = read_statistics('output-cutoff-c1/statistics')
cutoffc1_oscillating_data = [
    'Cutoff-c1, oscillating V',
    cutoffc1_df['Cell Score Standard Deviation: '].mean(),
    cutoffc1_df['Average particle distribution score: '].mean(),
    cutoffc1_df['Maximal particle distribution score: '].mean()
]

# Uniform kernel function, oscillating velocity
uniform_df = read_statistics('output-uniform/statistics')
uniform_oscillating_data = [
    'Uniform, oscillating V',
    uniform_df['Cell Score Standard Deviation: '].mean(),
    uniform_df['Average particle distribution score: '].mean(),
    uniform_df['Maximal particle distribution score: '].mean()
]

# Triangular kernel function, oscillating velocity
triangular_df = read_statistics('output-triangular/statistics')
triangular_oscillating_data = [
    'Triangular, oscillating V',
    triangular_df['Cell Score Standard Deviation: '].mean(),
    triangular_df['Average particle distribution score: '].mean(),
    triangular_df['Maximal particle distribution score: '].mean()
]

# Random kernel function, constant velocity
random_constant_df = read_statistics('output-random-constant-velocity/statistics')
random_constant_data = [
    'Random, constant V',
    random_constant_df['Cell Score Standard Deviation: '].mean(),
    random_constant_df['Average particle distribution score: '].mean(),
    random_constant_df['Maximal particle distribution score: '].mean()
]

# Gaussian kernel function, constant velocity
gaussian_constant_df = read_statistics('output-gaussian-constant-velocity/statistics')
gaussian_constant_data = [
    'Gaussian, constant V',
    gaussian_constant_df['Cell Score Standard Deviation: '].mean(),
    gaussian_constant_df['Average particle distribution score: '].mean(),
    gaussian_constant_df['Maximal particle distribution score: '].mean()
]

# Cutoff-w1 kernel function, constant velocity
cutoffw1_constant_df = read_statistics('output-cutoff-w1-constant-velocity/statistics')
cutoffw1_constant_data = [
    'Cutoff-w1, constant V',
    cutoffw1_constant_df['Cell Score Standard Deviation: '].mean(),
    cutoffw1_constant_df['Average particle distribution score: '].mean(),
    cutoffw1_constant_df['Maximal particle distribution score: '].mean(),
]

# Cutoff-c1 kernel function, constant velocity
cutoffc1_constant_df = read_statistics('output-cutoff-c1-constant-velocity/statistics')
cutoffc1_constant_data = [
    'Cutoff-c1, constant V',
    cutoffc1_constant_df['Cell Score Standard Deviation: '].mean(),
    cutoffc1_constant_df['Average particle distribution score: '].mean(),
    cutoffc1_constant_df['Maximal particle distribution score: '].mean()
]

# Uniform kernel function, constant velocity
uniform_constant_df = read_statistics('output-uniform-constant-velocity/statistics')
uniform_constant_data = [
    'Uniform, constant V',
    uniform_constant_df['Cell Score Standard Deviation: '].mean(),
    uniform_constant_df['Average particle distribution score: '].mean(),
    uniform_constant_df['Maximal particle distribution score: '].mean(),
]

# Triangular kernel function, constant velocity
triangular_constant_df = read_statistics('output-triangular-constant-velocity/statistics')
triangular_constant_data = [
    'Triangular, constant V',
    triangular_constant_df['Cell Score Standard Deviation: '].mean(),
    triangular_constant_df['Average particle distribution score: '].mean(),
    triangular_constant_df['Maximal particle distribution score: '].mean()
]

output_data_array = [
    random_oscillating_data,
    gaussian_oscillating_data,
    cutoffw1_oscillating_data,
    cutoffc1_oscillating_data,
    uniform_oscillating_data,
    triangular_oscillating_data,
    random_constant_data,
    gaussian_constant_data,
    cutoffw1_constant_data,
    cutoffc1_constant_data,
    uniform_constant_data,
    triangular_constant_data
]

column_names_output = [
    'Particle Removal Algorithm',
    'Time Averaged Cell Score Standard Deviation',
    'Time Averaged Mean Score',
    'Time Averaged Maximum Score',
]

output_dataframe = pd.DataFrame(output_data_array,columns=column_names_output)
output_dataframe.to_csv('removal_algorithm_comparison_data.csv',index=False)
