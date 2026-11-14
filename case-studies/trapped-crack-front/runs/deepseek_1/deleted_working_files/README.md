# Working files of deepseek_1 that the run deleted before it finished

The run wrote and ran these five scripts, then removed them with the command below. They are restored here from the
event stream of the run (the content it wrote, and what each printed when the run executed it), as a record of how the
run chose its response coefficients. `compute.py` and `scan.py` compare several sign patterns of the nonlocal term.

    cd /work && rm -f compute.py scan.py final.py check.py verify.py && python3 solve.py > solve_log.txt 2>&1 && python3 check_fixed.py && echo "--- files ---" && ls -la && echo "--- results.json ---" && cat results.json
