
# delete the generated polars
# it is useful if you want to benchmark the polar generation or if you have changed parameters

file1="data/ram_air_kite_foil_cl_polar.csv"
file2="data/ram_air_kite_foil_cd_polar.csv"
file3="data/ram_air_kite_foil_cm_polar.csv"
isfile(file1) && rm(file1)
isfile(file2) && rm(file2)
isfile(file3) && rm(file3)
