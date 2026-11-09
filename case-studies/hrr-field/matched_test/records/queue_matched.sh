#!/bin/bash
# Matched test queue (log P251): wave r1 (27 runs), then wave r2. 50 min per run.
H=<library-pilot-folder>/pilot_entry_format
export ROOT=<matched-runs>/runs LOGS=<matched-runs>-logs
for rep in r1 r2; do
  names=$(ls $ROOT | grep -- "-$rep-")
  for n in $names; do bash $H/launch_run.sh $n 3000; sleep 2; done
  until [ $(for n in $names; do [ -f $LOGS/$n.exit ] && echo x; done | wc -l) -ge $(echo $names | wc -w) ]; do sleep 30; done
  echo "wave $rep done $(date -u +%FT%TZ)"
done
