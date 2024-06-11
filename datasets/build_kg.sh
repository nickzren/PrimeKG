#!/bin/bash

# Log file
LOGFILE="build_kg.log"

# Recreate log file
echo "Starting KG pipeline" > $LOGFILE

# Function to run a script and check for errors
run_script() {
  local script_name=$1
  echo "Starting $script_name" | tee -a $LOGFILE
  python $script_name >> $LOGFILE 2>&1
  if [ $? -ne 0 ]; then
    echo "Error running $script_name" | tee -a $LOGFILE
    exit 1
  fi
}

# Run the sequence of Python scripts to create the final KG data
run_script "build_raw_kg.py"
run_script "build_idx2group.py"
run_script "build_embeds.py"
run_script "build_kg_grouped_diseases.py"
run_script "build_final_kg.py"

echo "KG pipeline completed successfully" | tee -a $LOGFILE