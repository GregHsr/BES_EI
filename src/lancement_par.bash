#!/bin/bash

make clean 
rm -f *.scl *.vec *.dat sortieEnsight.* script.log

# Définition des listes
list_re=(100 200 500 1000 5000 10000)
list_nx=(100 150 200 300)
list_schema=(1 2)

istock=10000
case_file="sortieEnsight.case"
template_file="data_template.txt"

export template_file istock case_file

# Fonction pour exécuter un calcul individuel
do_calculation() {
    re=$1
    nx=$2
    schema=$3

    output_file="data_${re}_${nx}_${schema}.txt"
    output_dir="Re_${re}_Nx_${nx}_Sc_${schema}"

    # Remplacement des valeurs et création du fichier de configuration
    sed -e "s/@Re@/$re/g" \
        -e "s/@Nx@/$nx/g" \
        -e "s/@Schema@/$schema/g" \
        "$template_file" > "$output_file"

    sed -e "s/@Re@/$re/g" \
        -e "s/@Nx@/$nx/g" \
        -e "s/@Schema@/$schema/g" \
        "$template_file" > "data.txt"

    echo "Fichier de configuration généré : $output_file" >> script.log

    # Création du dossier de sortie
    mkdir -p "$output_dir"

    cp *.f90 "$output_dir/"
    mv "$output_file" "$output_dir/"
    cp Makefile "$output_dir/"

    cd "$output_dir/"
    cp "$output_file" data.txt

    make

    # Exécution de cavite.exe
    start=$(date +%s)
    output=$(./cavite.exe)
    end=$(date +%s)
    runtime=$((end-start))

    echo "$output" >> script.log
    echo "Temps d'exécution : $runtime secondes" >> script.log

    iterations=$(echo "$output" | grep -oP 'Convergence atteinte en\s+\K[0-9]+')
    time_step_stockage=$(echo "$output" | grep -oP 'Time step stockage\s+\K[0-9]+')

    if [[ ! -f $case_file ]]; then
        echo "Le fichier $case_file n'existe pas." >> script.log
        exit 1
    fi

    number_of_steps=$((iterations / istock))
    last_time=$((number_of_steps * time_step_stockage))
    nb_st=$((number_of_steps + 1))
    last_line=$((18 + nb_st))

    sed -i "s/number of steps: *[0-9]\+/number of steps:     $number_of_steps/" "$case_file"
    sed -i "${last_line},\$d" "$case_file"

    cd ..
}

# Exporter la fonction pour que parallel puisse l'utiliser
export -f do_calculation

# Lancer les calculs en parallèle
echo "Lancement des calculs en parallèle..."
parallel do_calculation ::: ${list_re[@]} ::: ${list_nx[@]} ::: ${list_schema[@]}

make clean
rm data.txt
