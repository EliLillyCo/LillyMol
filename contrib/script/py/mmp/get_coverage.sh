#!/usr/bin/env bash
coverage erase
coverage run mmp_data_objects.py
coverage run -a mmp_dicer_functions.py
coverage run -a mmp_enum_mols_from_pairs.py
coverage run -a mmp_math_functions.py
coverage run -a mmp_mcss_objects.py
coverage run -a mmp_objects.py
coverage run -a mmp_pairs_objects.py
coverage run -a mmp_series_object.py
coverage run -a mmp_stats_functions.py
coverage run -a mmp_stats_functions_timer.py
coverage report -m mmp_data_objects.py mmp_dicer_functions.py mmp_enum_mols_from_pairs.py mmp_math_functions.py mmp_mcss_objects.py mmp_objects.py mmp_pairs_objects.py mmp_series_object.py mmp_stats_functions.py mmp_stats_functions_timer.py
