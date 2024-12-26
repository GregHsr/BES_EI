make clean 
rm -f *.scl *.vec *.dat sortieEnsight.* 
make

#!/bin/bash

# Définition des listes
list_re=(100 , 200)
list_nx=(100)
list_schema=(1)

istock=500
case_file="sortieEnsight.case"

# Nom du fichier modèle
template_file="data_template.txt"

# Boucle sur toutes les combinaisons
for re in "${list_re[@]}"; do
  for nx in "${list_nx[@]}"; do
    for schema in "${list_schema[@]}"; do
      
      # Nom du fichier de configuration et du dossier de sortie
      output_file="data_${re}_${nx}_${schema}.txt"
      output_dir="${re}_${nx}_${schema}"

      # Remplacement des valeurs et création du fichier de configuration
      sed -e "s/@Re@/$re/g" \
          -e "s/@Nx@/$nx/g" \
          -e "s/@Schema@/$schema/g" \
          "$template_file" > "$output_file"

      sed -e "s/@Re@/$re/g" \
          -e "s/@Nx@/$nx/g" \
          -e "s/@Schema@/$schema/g" \
          "$template_file" > "data.txt"

      echo "Fichier de configuration généré : $output_file"

      # Création du dossier de sortie
      mkdir -p "$output_dir"
      
      # Exécution de cavite.exe et capture de la sortie
      output=$(./cavite.exe)

      echo "$output"
      
      # Extraction des valeurs nécessaires de la sortie
      iterations=$(echo "$output" | grep -oP 'Convergence atteinte en\s+\K[0-9]+')
      time_step_stockage=$(echo "$output" | grep -oP 'Time step stockage\s+\K[0-9]+')

      if [[ ! -f $case_file ]]; then
        echo "Le fichier $case_file n'existe pas."
        exit 1
      fi
      
      # Calcul du nombre de steps et de l'incrément
      number_of_steps=$((iterations / istock))
      last_time=$((number_of_steps * time_step_stockage))
      nb_st=$(($number_of_steps + 1))
      last_line=$((18 + $nb_st))

      echo "Nombre d'itérations : $iterations"
      echo "Nombre de steps : $number_of_steps"
      echo "Dernier temps : $last_time"

      sed -i "s/number of steps: *[0-9]\+/number of steps:     $number_of_steps/" "$case_file"
      sed -i "${last_line},\$d" "$case_file"

      # Déplacement des fichiers générés dans le dossier de sortie
      for ext in "*.scl" "*.vec" "*.dat" "sortieEnsight.*" "data_*_*_*.txt"; do
        mv $ext "$output_dir/" 2>/dev/null || echo "Aucun fichier $ext à déplacer pour $output_file"
      done
      echo "Fichiers déplacés dans le dossier : $output_dir"
    done
  done
done

make clean
rm data.txt