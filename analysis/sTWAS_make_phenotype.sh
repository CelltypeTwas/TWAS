

mkdir -p /data && cd /data
mkdir -p /sh && cd /sh

for bamfiles in `ls /sudo_bam/*.bam`; do
    echo "regtools junctions extract -a 5 -m 30 -M 500000 -s XS $bamfiles -o ${bamfiles}.junc" >> all
done

mkdir -p split_sh
input_file="all"
output_dir="split_sh"
mkdir -p "$output_dir"
output_prefix="$output_dir/small_file_"
split -l 1111 --additional-suffix=.pbs "$input_file" "$output_prefix"






