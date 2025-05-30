#!/bin/bash

logfile="/Users/gc3045/scmail_v1/LAML/example2_bothscores.log"

# For llh diff
avg_llh_diff=$(grep "llh diff between EMsolvers (old - new):" "$logfile" | awk '{sum += $NF; count++} END {print sum/count}')

# For time diff
avg_time_diff=$(grep "time diff between EMsolvers (old - new):" "$logfile" | awk '{sum += $NF; count++} END {print sum/count}')

# For param (phi) diff
avg_phi_diff=$(grep "param (phi) diff between EMsolvers (old - new):" "$logfile" | awk '{sum += $NF; count++} END {print sum/count}')

# For param (nu) diff
avg_nu_diff=$(grep "param (nu) diff between EMsolvers (old - new):" "$logfile" | awk '{sum += $NF; count++} END {print sum/count}')

# Print results
echo "Average llh diff: $avg_llh_diff"
echo "Average time diff: $avg_time_diff"
echo "Average phi diff: $avg_phi_diff"
echo "Average nu diff: $avg_nu_diff"

