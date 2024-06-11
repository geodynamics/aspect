import matplotlib.pyplot as plt
import numpy as np
import math

# Parametri
r_inner = 3481000  # Raggio interno in metri
r_outer = 6371000  # Raggio esterno in metri
p_large = 700.0  # Ampiezza delle grandi perturbazioni
p_small = 30.0  # Ampiezza delle piccole perturbazioni
frequency_large = 3  # Frequenza delle grandi perturbazioni
frequency_small = 30  # Frequenza delle piccole perturbazioni
spread_angle_large_degrees = 90  # Angolo di diffusione delle grandi perturbazioni in gradi
start_angle_degrees = 270  # Angolo di partenza per la perturbazione in gradi

# Funzione di temperatura modificata per utilizzare i gradi e gestire gli array
def temperature_function_degrees_array(r, theta_degrees):
    # Aggiusta theta per spostare l'angolo di inizio della perturbazione
    adjusted_theta_degrees = theta_degrees + start_angle_degrees
    adjusted_theta_degrees = adjusted_theta_degrees % 360

    # Inizializza l'array della temperatura
    temp = np.zeros_like(r)

    # Determina la grande perturbazione usando un approccio a funzione gradino
    large_perturbation_mask = (adjusted_theta_degrees < spread_angle_large_degrees) | (adjusted_theta_degrees > 360 - spread_angle_large_degrees)
    temp[large_perturbation_mask] = p_large * np.sin(np.pi * (r[large_perturbation_mask] - r_inner) / (r_outer - r_inner) + math.pi) * np.sin(frequency_large * np.radians(theta_degrees[large_perturbation_mask]))

    # Applica piccole perturbazioni dove non si applicano le grandi perturbazioni
    small_perturbation_mask = ~large_perturbation_mask
    temp[small_perturbation_mask] = p_small * np.sin(np.pi * (r[small_perturbation_mask] - r_inner) / (r_outer - r_inner)) * np.sin(frequency_small * np.radians(theta_degrees[small_perturbation_mask]))

    return temp
    
    
    
    
# Generazione dei dati per la visualizzazione
theta = np.linspace(0, 360, 360)  # Angoli da 0 a 360 gradi
radius = np.linspace(r_inner, r_outer, 100)  # Raggi da r_inner a r_outer

# Creazione di una griglia per il plotting
R, Theta = np.meshgrid(radius, theta)
Temp = temperature_function_degrees_array(R, Theta)

def write_ascii_data(radius, theta_degrees, temperatures, filename):
    with open(filename, 'w') as file:
        file.write("# Test data for ascii data initial conditions.\n")
        file.write("# Only next line is parsed in format: [nx] [ny] [nz] because of keyword \"POINTS:\"\n")
        file.write("# POINTS: 100 360 1\n")
        file.write("# Columns: r phi      temperature [K]\n")

        for j in range(len(theta_degrees)):
            for i in range(len(radius)):
                r = radius[i]
                phi = theta_degrees[j]
                temp = temperatures[j, i]

                # Formattazione migliorata
                r_str = f"{r:.1f}"
                phi_str = f"{phi:.4f}" if not phi.is_integer() else f"{int(phi)}"
                temp_str = f"{temp:.1f}" if not temp.is_integer() else f"{int(temp)}"

                # Allineamento delle colonne con spostamento della seconda colonna
                file.write(f"{r_str:<10} {phi_str:<11}{temp_str}\n")


# Conversione delle coordinate polari in coordinate cartesiane per il plotting
X = R * np.cos(np.radians(Theta))
Y = R * np.sin(np.radians(Theta))

# Plotting
plt.figure(figsize=(8, 6))
contour = plt.contourf(X, Y, Temp, levels=50, cmap='viridis')
plt.colorbar(contour, label='Temperatura')
plt.title('Distribuzione della Temperatura in Coordinate Polari')
plt.xlabel('Coordinata X (m)')
plt.ylabel('Coordinata Y (m)')
plt.axis('equal')
plt.show()


# Uso della funzione write_ascii_data per salvare i dati
write_ascii_data(radius, theta, Temp, 'perturbation_ascii.txt')

import math

def convert_to_radians(degrees):
    """Convert degrees to radians."""
    return degrees * (math.pi / 180)

def convert_file_angles(input_file_path, output_file_path):
    """Convert the angles in the second column of the file from degrees to radians."""
    with open(input_file_path, 'r') as file:
        content = file.readlines()

    converted_lines = content[:4]  # Keep the header lines unchanged
    for line in content[4:]:
        columns = line.split()
        if len(columns) == 3:
            radius, angle_deg, temperature = columns
            angle_rad = convert_to_radians(float(angle_deg))
            converted_line = f"{radius}  {angle_rad:.8f}  {temperature}\n"
            converted_lines.append(converted_line)
        else:
            converted_lines.append(line)

    with open(output_file_path, 'w') as file:
        file.writelines(converted_lines)

# Example usage
input_file_path = 'perturbation_ascii.txt'  # Replace with your input file path
output_file_path = 'perturbation_ascii.txt'  # Replace with your desired output file path

convert_file_angles(input_file_path, output_file_path)
