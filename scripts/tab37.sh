#!/bin/bash

CC="clang"

folder="/Users/antonmosunov/Desktop/test"
# Make sure that discriminant_bound is divisible by 15
discriminant_bound=$((15 * 2**37))
number_of_files=$((2**12))
# Bitsize values in order: 23, 47, 95 mod 120,   7, 15 mod 24,   4, 8 mod 16,   3 mod 8
bitsize_values=(25 25 25 26 25 25 24 27)
degree_after_bundling=$((2**26))
degree_values=(
    $((discriminant_bound / 120))
    $((discriminant_bound / 120))
    $((discriminant_bound / 120))
    $((discriminant_bound / 24))
    $((discriminant_bound / 24))
    $((discriminant_bound / 16))
    $((discriminant_bound / 16))
    $((discriminant_bound / 8))
)
bundle_values=(
    $((degree_values[0] / degree_after_bundling))
    $((degree_values[1] / degree_after_bundling))
    $((degree_values[2] / degree_after_bundling))
    $((degree_values[3] / degree_after_bundling))
    $((degree_values[4] / degree_after_bundling))
    $((degree_values[5] / degree_after_bundling))
    $((degree_values[6] / degree_after_bundling))
    $((degree_values[7] / degree_after_bundling))
)

set -x

# echo "${degree_values[@]}"
# echo "${bundle_values[@]}"
# echo "${bitsize_values[@]}"

# exit 0

#0. Compiling all files
cd ..
$CC -fopenmp -lflint *.c -L/usr/local/lib -I/usr/local/include -o polymult
cd post-processing
$CC -fopenmp -DDELETE_INPUT_FILES sum.c -o sum
$CC -fopenmp -DDELETE_INPUT_FILES merge.c -o merge
$CC -fopenmp -DDELETE_INPUT_FILES scale.c -o scale
$CC replace.c -o replace
cd ..





#1. Generating the polynomial sum H(120k+23)q^k 
date
polynomial_degree=${degree_values[0]}
bundle=${bundle_values[0]}
bitsize=${bitsize_values[0]}

    #1.1 Generating the first summand
summand1_name="h23mod120summand1."
./polymult $polynomial_degree $number_of_files $bundle $bitsize $summand1_name $folder   1 0 1 1 3   2 0 2 2 15   1 0 2 1 3   2 1 2 8 15   1 0 2 1 3   2 1 2 7 15   1 0 2 2 3   2 3 2 13 15   1 0 2 2 3   2 3 6 4 5   1 0 6 1 1   1 0 6 1 5   1 0 6 0 1

    #1.2 Generating the second summand
summand2_name="h23mod120summand2."
./polymult $polynomial_degree $number_of_files $bundle $bitsize $summand2_name $folder   1 0 3 1 1   2 1 2 2 15   1 0 6 1 1   2 2 2 8 15   1 0 6 1 1   1 1 2 7 15   1 0 6 0 1   1 3 2 13 15   1 0 6 0 1

    #1.3 Adding the two polynomials
result_name="h23mod120."
output="$folder/$result_name"
input1="$folder/$summand1_name"
input2="$folder/$summand2_name"
./post-processing/sum $number_of_files $output $input1 $input2

    #1.4 Deleting temporary files
# rm "$folder/$summand1_name"*
# rm "$folder/$summand2_name"*





#2. Generating the polynomial sum H(120k+47)q^k
date
polynomial_degree=${degree_values[1]}
bundle=${bundle_values[1]}
bitsize=${bitsize_values[1]}

    #2.1 Generating the first summand
summand1_name="h47mod120summand1."
./polymult $polynomial_degree $number_of_files $bundle $bitsize $summand1_name $folder   1 0 1 1 3   2 0 2 4 15   1 0 2 1 3   2 3 2 14 15   1 0 2 1 3   2 0 2 1 15   1 0 2 2 3   2 2 2 11 15   1 0 2 2 3   2 1 6 2 5   1 0 6 1 1   1 1 6 3 5   1 0 6 0 1

    #2.2 Generating the second summand
summand2_name="h47mod120summand2."
./polymult $polynomial_degree $number_of_files $bundle $bitsize $summand2_name $folder   1 0 3 1 1   2 1 2 4 15   1 0 6 1 1   2 4 2 14 15   1 0 6 1 1   1 0 2 1 15   1 0 6 0 1   1 2 2 11 15   1 0 6 0 1

    #2.3 Adding the two polynomials
result_name="h47mod120."
output="$folder/$result_name"
input1="$folder/$summand1_name"
input2="$folder/$summand2_name"
./post-processing/sum $number_of_files $output $input1 $input2

    #2.4 Deleting temporary files
# rm "$folder/$summand1_name"*
# rm "$folder/$summand2_name"*





#3. Generating the polynomial sum H(120k+95)q^k
date
polynomial_degree=${degree_values[2]}
bundle=${bundle_values[2]}
bitsize=${bitsize_values[2]}
factor=1

if (($factor == 1)); then
    bitsize=$((bitsize - 1))
fi

    #3.1 Generating the first summand
summand1_name="h95mod120summand1."
./polymult $polynomial_degree $number_of_files $bundle $bitsize $summand1_name $folder   $factor 0 1 1 3   2 1 10 2 3   1 0 2 1 3   2 0 10 1 3   1 0 2 2 3   1 0 30 0 1   1 0 6 1 1   1 3 30 1 1   1 0 6 0 1

    #3.2 Generating the second summand
summand2_name="h95mod120summand2."
./polymult $polynomial_degree $number_of_files $bundle $bitsize $summand2_name $folder   $factor 0 3 1 1   2 2 10 2 3   1 0 6 1 1   1 0 10 1 3   1 0 6 0 1

    #3.3 Adding the two polynomials
result_name="h95mod120."
output="$folder/$result_name"
input1="$folder/$summand1_name"
input2="$folder/$summand2_name"
./post-processing/sum $number_of_files $output $input1 $input2

if (($factor == 1)); then
    ./post-processing/scale $number_of_files $output $output multiply 2
fi

    #3.4 Deleting temporary files
# rm "$folder/$summand1_name"*
# rm "$folder/$summand2_name"*





#4. Producing the file containing H(24k+23) (0s are recorded for 71 (mod 120) and 119 (mod 120)) and deleting all files containing mod 120 data
date

./post-processing/merge $folder $number_of_files 23 24 120
# result_name="h23mod120."
# rm "$folder/$result_name"*
# result_name="h47mod120."
# rm "$folder/$result_name"*
# result_name="h95mod120."
# rm "$folder/$result_name"*





#5. Generating the polynomial sum H(24k+7)q^k
date
polynomial_degree=${degree_values[3]}
bundle=${bundle_values[3]}
bitsize=${bitsize_values[3]}

    #5.1 Generating the polynomial 2*sum H(24k+7)q^k
result_name="h7mod24."
./polymult $polynomial_degree $number_of_files $bundle $bitsize $result_name $folder   1 0 1 1 1   1 0 1 1 3   1 0 1 1 1   1 0 4 1 3   1 0 4 0 1   2 1 4 2 3   1 0 4 1 1

    #5.2 Dividing by 2
output="$folder/$result_name"
input="$folder/$result_name"
./post-processing/scale $number_of_files $output $input divide 2





#6. Generating the polynomial sum H(24k+15)q^k
date
polynomial_degree=${degree_values[4]}
bundle=${bundle_values[4]}
bitsize=${bitsize_values[4]}

result_name="h15mod24."
./polymult $polynomial_degree $number_of_files $bundle $bitsize $result_name $folder   1 0 1 1 1   1 0 3 1 1   1 0 1 1 1   1 1 12 1 1   1 0 4 0 1   1 0 12 0 1   1 0 4 1 1





#7. Producing the file containing H(8k+7) (0s are recorded for 71 (mod 120) and 119 (mod 120)) and deleting all files containing mod 24 data
date

./post-processing/merge $folder $number_of_files 7 8 24
# result_name="h7mod24."
# rm "$folder/$result_name"*
# result_name="h15mod24."
# rm "$folder/$result_name"*
# result_name="h23mod24."
# rm "$folder/$result_name"*





#8. Generating the polynomial sum H(16k+4)q^k
date
polynomial_degree=${degree_values[5]}
bundle=${bundle_values[5]}
bitsize=${bitsize_values[5]}

    #8.1 Generating the polynomial 2*sum H(16k+4)q^k
result_name="h4mod16."
./polymult $polynomial_degree $number_of_files $bundle $bitsize $result_name $folder 1 0 2 1 1   1 0 2 0 1   1 0 2 0 1

    #8.2 Dividing by 2
output="$folder/$result_name"
input="$folder/$result_name"
./post-processing/scale $number_of_files $output $input divide 2

    #8.3 Replacing the 0th entry in h4mod16.0 with 1 (because Q(sqrt(-4)) has non-trivial automorphisms)
input="$folder/$result_name"
input+="0"
./post-processing/replace $input 0 1 





#9. Generating the polynomial sum H(16k+8)q^k
date
polynomial_degree=${degree_values[6]}
bundle=${bundle_values[6]}
bitsize=${bitsize_values[6]}

result_name="h8mod16."
./polymult $polynomial_degree $number_of_files $bundle $bitsize $result_name $folder 1 0 2 0 1   1 0 2 1 1   1 0 2 1 1





#10. Generating the polynomial sum H(8k+3)q^k
date
polynomial_degree=${degree_values[7]}
bundle=${bundle_values[7]}
bitsize=${bitsize_values[7]}

    #10.1 Generating the polynomial 3*sum H(8k+3)q^k
result_name="h3mod8."
./polymult $polynomial_degree $number_of_files $bundle $bitsize $result_name $folder 1 0 1 1 1   1 0 1 1 1   1 0 1 1 1

    #10.2 Dividing by 3
output="$folder/$result_name"
input="$folder/$result_name"
./post-processing/scale $number_of_files $output $input divide 3

    #10.3 Replacing the 0th entry in h3mod8.0 with 1 (because Q(sqrt(-3)) has non-trivial automorphisms)
input="$folder/$result_name"
input+="0"
./post-processing/replace $input 0 1 

date
