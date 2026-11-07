# REPRODUCE.md — reviewer checks

Single command (this folder):

    timeout 60 python3 review_check_independent.py

Writes `review_check_independent_output.txt` and `review_check_independent_results.json`.
Does not import `submission/` or `reference/` solvers. Reads those JSON files only for numeric comparison.
