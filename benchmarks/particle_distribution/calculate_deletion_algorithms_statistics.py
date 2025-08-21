import pandas as pd

# There may be a way to read the column names from the statistics file, but this will work unless the postprocessors are changed in the future.
column_names = [
    'Time step number',
    'Time (seconds)',
    'Time step size (seconds)',
    'Number of mesh cells',
    'Number of Stokes degrees of freedom',
    'Number of temperature degrees of freedom',
    'Number of degrees of freedom for all compositions',
    'Iterations for temperature solver',
    'Iterations for composition solver 1',
    'RMS velocity (m/s)',   
    'Max. velocity (m/s)',
    'Minimal temperature (K)', 
    'Average temperature (K)',
    'Maximal temperature (K)', 
    'Outward heat flux through boundary with indicator 0 ("left") (W)',
    'Outward heat flux through boundary with indicator 1 ("right") (W)', 
    'Outward heat flux through boundary with indicator 2 ("bottom") (W)',
    'Outward heat flux through boundary with indicator 3 ("top") (W)', 
    'Visualization file name',
    'Number of advected particles', 
    'Particle file name',
    'Minimal particle distribution score', 
    'Average particle distribution score',
    'Maximal particle distribution score', 
    'Minimum PDF standard deviation',
    'Mean of PDF standard deviation', 
    'Maximum PDF standard deviation',
    'Mean of PDF minimum values',  
    'Mean PDF maximum values',
    'Minimum of PDF minimum values',
    'Maximum of PDF maximum values'
]

# Random selection / no kernel function, oscillating velocity
random_df = pd.read_table('output-random/statistics',names=column_names,comment='#',engine='python',sep='\s+')
random_oscillating_data = [
    'Random, oscillating V',
    random_df['Minimal particle distribution score'].mean(),
    random_df['Average particle distribution score'].mean(),
    random_df['Maximal particle distribution score'].mean(),
    random_df['Minimum PDF standard deviation'].mean(),
    random_df['Mean of PDF standard deviation'].mean(),
    random_df['Maximum PDF standard deviation'].mean(),
    random_df['Mean of PDF minimum values'].mean(),
    random_df['Mean PDF maximum values'].mean(),
    random_df['Minimum of PDF minimum values'].min(),
    random_df['Maximum of PDF maximum values'].max()
]

# Gaussian kernel function, oscillating velocity
gaussian_df = pd.read_table('output-gaussian/statistics',names=column_names,comment='#',engine='python',sep='\s+')
gaussian_oscillating_data = [
    'Gaussian, oscillating V',
    gaussian_df['Minimal particle distribution score'].mean(),
    gaussian_df['Average particle distribution score'].mean(),
    gaussian_df['Maximal particle distribution score'].mean(),
    gaussian_df['Minimum PDF standard deviation'].mean(),
    gaussian_df['Mean of PDF standard deviation'].mean(),
    gaussian_df['Maximum PDF standard deviation'].mean(),
    gaussian_df['Mean of PDF minimum values'].mean(),
    gaussian_df['Mean PDF maximum values'].mean(),
    gaussian_df['Minimum of PDF minimum values'].mean(),
    gaussian_df['Maximum of PDF maximum values'].mean()
]

# Cutoff-w1 kernel function, oscillating velocity
cutoffw1_df = pd.read_table('output-cutoff-w1/statistics',names=column_names,comment='#',engine='python',sep='\s+')
cutoffw1_oscillating_data = [
    'Cutoff-w1, oscillating V',
    cutoffw1_df['Minimal particle distribution score'].mean(),
    cutoffw1_df['Average particle distribution score'].mean(),
    cutoffw1_df['Maximal particle distribution score'].mean(),
    cutoffw1_df['Minimum PDF standard deviation'].mean(),
    cutoffw1_df['Mean of PDF standard deviation'].mean(),
    cutoffw1_df['Maximum PDF standard deviation'].mean(),
    cutoffw1_df['Mean of PDF minimum values'].mean(),
    cutoffw1_df['Mean PDF maximum values'].mean(),
    cutoffw1_df['Minimum of PDF minimum values'].mean(),
    cutoffw1_df['Maximum of PDF maximum values'].mean()
]

# Uniform kernel function, oscillating velocity
uniform_df = pd.read_table('output-uniform/statistics',names=column_names,comment='#',engine='python',sep='\s+')
uniform_oscillating_data = [
    'Uniform, oscillating V',
    uniform_df['Minimal particle distribution score'].mean(),
    uniform_df['Average particle distribution score'].mean(),
    uniform_df['Maximal particle distribution score'].mean(),
    uniform_df['Minimum PDF standard deviation'].mean(),
    uniform_df['Mean of PDF standard deviation'].mean(),
    uniform_df['Maximum PDF standard deviation'].mean(),
    uniform_df['Mean of PDF minimum values'].mean(),
    uniform_df['Mean PDF maximum values'].mean(),
    uniform_df['Minimum of PDF minimum values'].mean(),
    uniform_df['Maximum of PDF maximum values'].mean()
]

# Triangular kernel function, oscillating velocity
triangular_df = pd.read_table('output-triangular/statistics',names=column_names,comment='#',engine='python',sep='\s+')
triangular_oscillating_data = [
    'Triangular, oscillating V',
    triangular_df['Minimal particle distribution score'].mean(),
    triangular_df['Average particle distribution score'].mean(),
    triangular_df['Maximal particle distribution score'].mean(),
    triangular_df['Minimum PDF standard deviation'].mean(),
    triangular_df['Mean of PDF standard deviation'].mean(),
    triangular_df['Maximum PDF standard deviation'].mean(),
    triangular_df['Mean of PDF minimum values'].mean(),
    triangular_df['Mean PDF maximum values'].mean(),
    triangular_df['Minimum of PDF minimum values'].mean(),
    triangular_df['Maximum of PDF maximum values'].mean()
]

# Random kernel function, constant velocity
random_constant_df = pd.read_table('output-random-constant-velocity/statistics',names=column_names,comment='#',engine='python',sep='\s+')
random_constant_data = [
    'Random, constant V',
    random_constant_df['Minimal particle distribution score'].mean(),
    random_constant_df['Average particle distribution score'].mean(),
    random_constant_df['Maximal particle distribution score'].mean(),
    random_constant_df['Minimum PDF standard deviation'].mean(),
    random_constant_df['Mean of PDF standard deviation'].mean(),
    random_constant_df['Maximum PDF standard deviation'].mean(),
    random_constant_df['Mean of PDF minimum values'].mean(),
    random_constant_df['Mean PDF maximum values'].mean(),
    random_constant_df['Minimum of PDF minimum values'].min(),
    random_constant_df['Maximum of PDF maximum values'].max()
]

# Gaussian kernel function, constant velocity
gaussian_constant_df = pd.read_table('output-gaussian-constant-velocity/statistics',names=column_names,comment='#',engine='python',sep='\s+')
gaussian_constant_data = [
    'Gaussian, constant V',
    gaussian_constant_df['Minimal particle distribution score'].mean(),
    gaussian_constant_df['Average particle distribution score'].mean(),
    gaussian_constant_df['Maximal particle distribution score'].mean(),
    gaussian_constant_df['Minimum PDF standard deviation'].mean(),
    gaussian_constant_df['Mean of PDF standard deviation'].mean(),
    gaussian_constant_df['Maximum PDF standard deviation'].mean(),
    gaussian_constant_df['Mean of PDF minimum values'].mean(),
    gaussian_constant_df['Mean PDF maximum values'].mean(),
    gaussian_constant_df['Minimum of PDF minimum values'].min(),
    gaussian_constant_df['Maximum of PDF maximum values'].max()
]

# Cutoff-w1 kernel function, constant velocity
cutoffw1_constant_df = pd.read_table('output-cutoff-w1-constant-velocity/statistics',names=column_names,comment='#',engine='python',sep='\s+')
cutoffw1_constant_data = [
    'Cutoff-w1, constant V',
    cutoffw1_constant_df['Minimal particle distribution score'].mean(),
    cutoffw1_constant_df['Average particle distribution score'].mean(),
    cutoffw1_constant_df['Maximal particle distribution score'].mean(),
    cutoffw1_constant_df['Minimum PDF standard deviation'].mean(),
    cutoffw1_constant_df['Mean of PDF standard deviation'].mean(),
    cutoffw1_constant_df['Maximum PDF standard deviation'].mean(),
    cutoffw1_constant_df['Mean of PDF minimum values'].mean(),
    cutoffw1_constant_df['Mean PDF maximum values'].mean(),
    cutoffw1_constant_df['Minimum of PDF minimum values'].min(),
    cutoffw1_constant_df['Maximum of PDF maximum values'].max()
]

# Uniform kernel function, constant velocity
uniform_constant_df = pd.read_table('output-uniform-constant-velocity/statistics',names=column_names,comment='#',engine='python',sep='\s+')
uniform_constant_data = [
    'Uniform, constant V',
    uniform_constant_df['Minimal particle distribution score'].mean(),
    uniform_constant_df['Average particle distribution score'].mean(),
    uniform_constant_df['Maximal particle distribution score'].mean(),
    uniform_constant_df['Minimum PDF standard deviation'].mean(),
    uniform_constant_df['Mean of PDF standard deviation'].mean(),
    uniform_constant_df['Maximum PDF standard deviation'].mean(),
    uniform_constant_df['Mean of PDF minimum values'].mean(),
    uniform_constant_df['Mean PDF maximum values'].mean(),
    uniform_constant_df['Minimum of PDF minimum values'].min(),
    uniform_constant_df['Maximum of PDF maximum values'].max()
]

# Triangular kernel function, constant velocity
triangular_constant_df = pd.read_table('output-triangular-constant-velocity/statistics',names=column_names,comment='#',engine='python',sep='\s+')
triangular_constant_data = [
    'Triangular, constant V',
    triangular_constant_df['Minimal particle distribution score'].mean(),
    triangular_constant_df['Average particle distribution score'].mean(),
    triangular_constant_df['Maximal particle distribution score'].mean(),
    triangular_constant_df['Minimum PDF standard deviation'].mean(),
    triangular_constant_df['Mean of PDF standard deviation'].mean(),
    triangular_constant_df['Maximum PDF standard deviation'].mean(),
    triangular_constant_df['Mean of PDF minimum values'].mean(),
    triangular_constant_df['Mean PDF maximum values'].mean(),
    triangular_constant_df['Minimum of PDF minimum values'].min(),
    triangular_constant_df['Maximum of PDF maximum values'].max()
]

output_data_array = [
    random_oscillating_data,
    gaussian_oscillating_data,
    cutoffw1_oscillating_data,
    uniform_oscillating_data,
    triangular_oscillating_data,
    random_constant_data,
    gaussian_constant_data,
    cutoffw1_constant_data,
    uniform_constant_data,
    triangular_constant_data
]

column_names_output = [
    'Deletion Algorithm',
    'Minimum Score Mean',
    'Average Score Mean',
    'Maximum Score Mean',
    'Minimum Standard Deviation',
    'Average Standard Deviation',
    'Maximum Standard Deviation',
    'Average PDF Minimum',
    'Average PDF Maximum',
    'Minimum PDF Minimum',
    'Maximum PDF Maximum'
]

output_dataframe = pd.DataFrame(output_data_array,columns=column_names_output)
output_dataframe.to_csv('deletionalgorithm_comparison_data.csv',index=False)
