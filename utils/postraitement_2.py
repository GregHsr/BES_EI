#%% Import libraries
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

#%% Read data
# Read documentation values

u_ghia = pd.read_csv('Prof_u_Ghia_Re_100_a_10000.dat', sep='\s+')
y_ghia = u_ghia['y']

v_ghia = pd.read_csv('Prof_v_Ghia_Re_100_a_10000.dat', sep='\s+')
x_ghia = v_ghia['x']
print(v_ghia)

# Read simulation values
Re = [100, 400, 1000, 3200, 5000]
# Initialisation d'une liste pour stocker les colonnes correspondant à chaque Re
columns = {'x': None}  # On commence par une colonne 'x'

# Lecture des fichiers et construction des colonnes
for re in Re:
    # Chargement des données simulées
    u_simu = pd.read_csv(f'../src/Re_{re}_Nx_100_Sc_1/u.dat', sep='\s+')
    
    # Vérification que 'y' et 'u' existent dans le fichier
    if 'y' not in u_simu.columns or 'u' not in u_simu.columns:
        raise ValueError(f"Le fichier ../src/u_Re_{re}.dat doit contenir les colonnes 'y' et 'u'.")
    
    # La colonne 'y' sera utilisée comme 'x' si elle n'est pas encore définie
    if columns['x'] is None:
        columns['x'] = u_simu['y']
    
    # Ajout des données pour ce Re
    columns[f'Re={re}'] = u_simu['u']

# Création de la dataframe finale
frame_simu = pd.DataFrame(columns)

# Affichage pour vérification
print(frame_simu)

#%% Comparaison
# Liste des couleurs pour chaque Reynolds
colors = {
    100: 'blue',
    400: 'green',
    1000: 'orange',
    3200: 'red',
    5000: 'purple',
    7500: 'brown',
    10000: 'cyan'
}

# Conversion des colonnes en float (si nécessaire)
df_articles = u_ghia.astype(float)  # Supposons que df_articles contient vos données d'article
df_simu = frame_simu.astype(float)  # frame_simu contient vos données simulées

# Création du graphiqueimport matplotlib.pyplot as plt
plt.figure(figsize=(12, 8))

# Tracer les données des articles
for col in df_articles.columns[1:]:  # On ignore la colonne 'x'
    re_value = int(col.split('=')[1])  # Extraction de la valeur de Re
    if re_value in colors and re_value in Re:
        plt.plot(
            df_articles['y'], df_articles[col],
            label=f'Article Re={re_value}', color=colors[re_value], linestyle='-'
        )

# Tracer les données de simulation
for col in df_simu.columns[1:]:  # On ignore la colonne 'x'
    re_value = int(col.split('=')[1])  # Extraction de la valeur de Re
    if re_value in colors:
        plt.plot(
            df_simu['x'][::5], df_simu[col][::5],
            label=f'Simulation Re={re_value}', color=colors[re_value], linestyle='None', marker='x'
        )

# Configurations du graphique
plt.xlabel('x')
plt.ylabel('u')
plt.title('Comparaison des Données Simulées et Articles par Reynolds')
plt.legend()
plt.grid(True)
plt.show()


#u_simu = pd.read_csv('../src/u.dat', sep='\s+')
#y_simu = u_simu['y']
#v_simu = pd.read_csv('../src/v.dat', sep='\s+')
#x_simu = v_simu['x']

# #%% Analyse for Re=100
# 
# u_ghia_Re100 = u_ghia['Re=100']
# u_simu_Re100 = u_simu['u']
# v_ghia_Re100 = v_ghia['Re=100']
# v_simu_Re100 = v_simu['v']
# 
# 
# #%% Plot u velocity
# plt.figure()
# plt.plot(y_ghia, u_ghia_Re100, label='Ghia et al.')
# plt.plot(y_simu, u_simu_Re100, label='Simulation')
# plt.xlabel('y')
# plt.ylabel('u')
# plt.legend()
# plt.title('u velocity for Re=100')
# plt.show()
# 
# 
# #%% Plot v velocity
# plt.figure()
# plt.plot(x_ghia, v_ghia_Re100, label='Ghia et al.')
# plt.plot(x_simu, v_simu_Re100, label='Simulation')
# plt.xlabel('x')
# plt.ylabel('v')
# plt.legend()
# plt.title('v velocity for Re=100')
# plt.show()
