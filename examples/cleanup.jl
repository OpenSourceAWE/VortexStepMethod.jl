
# delete the generated polars
# this might be needed when the software or Julia have been updated
# it is also useful if you want to benchmark the polar generation

file1="data/ram_air_kite_body_info.bin"
file2="data/ram_air_kite_foil_polar.bin"
if isfile(file1)
    rm(file1)
end
if isfile(file2)
    rm(file2)
end