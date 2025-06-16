#!/bin/bash

# Run GBM43 or GBM6 analysis

# Path to the R script
# R_SCRIPT_PATH="./R/xenografts/xenograft_pipeline.R"
#
# cat "$R_SCRIPT_PATH" > gbm43_run_code.out
#
# # Run the Rscript in the background safely using nohup
# nohup Rscript "$R_SCRIPT_PATH" gbm43_1 TRUE 15 TRUE > gbm43_out_1.out &
# nohup Rscript "$R_SCRIPT_PATH" gbm43_2 FALSE 15 TRUE > gbm43_out_2.out &
# nohup Rscript "$R_SCRIPT_PATH" gbm43_3 TRUE NULL TRUE > gbm43_out_3.out &
# nohup Rscript "$R_SCRIPT_PATH" gbm43_4 FALSE NULL TRUE > gbm43_out_4.out &
# nohup Rscript "$R_SCRIPT_PATH" gbm43_5 TRUE 15 FALSE > gbm43_out_5.out &
# nohup Rscript "$R_SCRIPT_PATH" gbm43_6 FALSE 15 FALSE > gbm43_out_6.out &
# nohup Rscript "$R_SCRIPT_PATH" gbm43_7 TRUE NULL FALSE > gbm43_out_7.out &
# nohup Rscript "$R_SCRIPT_PATH" gbm43_8 FALSE NULL FALSE > gbm43_out_8.out &

# Run SB28 analysis
R_SCRIPT_PATH="../R/xenografts/xenograft_sb28_pipeline.R"

cat "$R_SCRIPT_PATH" > sb28_run_code.out

# Run the Rscript in the background safely using nohup
nohup Rscript "$R_SCRIPT_PATH" sb28_1 TRUE 15 FALSE > sb28_out_1.out &
