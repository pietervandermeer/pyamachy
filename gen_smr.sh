python generate_smr.py db -cnone --orbit=$1 $2_raw.h5
python generate_smr.py db -c1,2 --orbit=$1 $2_12.h5
python generate_smr.py db -c1,2,3 --orbit=$1 $2_123.h5

